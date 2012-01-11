/* Split scans of listing files into character cells.
 * Duncan Smith
 * 2009, 2011
 *
 * gcc ess-ocr.c -ggdb3 -Wall -lnetpbm -lm -std=c99
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include <pgm.h>

#define NO_ROW -1

// Looking for a division between lines
#define HORZ 0

// Looking for a division between characters
#define VERT 1

int *rowdown	(gray **image, gray max, int direction, int num_rows, int num_cols, float distance);

float *countup	(gray **image, gray max, int direction, int num_rows, int num_cols);
float row_count	(gray **image, gray max, int direction, int num_rows, int num_cols, int offset);

int print_char	(gray **image, gray max, int top, int bottom, int left, int right);
int query_char	(gray **image, gray max, int top, int bottom, int left, int right);
struct recognition *recognize_char	(gray **image, gray max, int top, int bottom, int left, int right);
int train_char	(gray **image, gray max, int top, int bottom, int left, int right, int codepoint);


#define MAX_RECOGNITION 16
struct recognition {
	struct recognized_char {
		int codepoint;
		float certainty;
	} guess[MAX_RECOGNITION];
};

struct trained_picture {
	gray **data;
	int example_count;
	struct trained_picture *next;
};

#define UNDEFINED_CODEPOINT 0
#define MAX_CODEPOINT MAXINT
struct trained_char {
	int codepoint;
	struct trained_picture *pictures;
};

struct trained_char *training_data;
int training_size = 0;

int main(int argc, char *argv[])
{
	gray **image;
	FILE *fp;
	gray max;
	int prows, pcols;
	int *row_loc;
	int *col_loc;

	int c_row, c_col;

	pgm_init(&argc, argv);

	fp = fopen(argv[1], "r");

	image = pgm_readpgm(fp, &pcols, &prows, &max);

	puts("Finding columns");
	col_loc = rowdown(image, max, VERT, prows, pcols, 12.75);
	puts("Finding rows");
	row_loc = rowdown(image, max, HORZ, prows, pcols, 21.5);

	printf("%d rows and %d columns\n", row_loc[0], col_loc[0]);

/*	for (c_row = 1; c_row < row_loc[0] - 1; c_row++) {
		for (c_col = 1; c_col < col_loc[0]; c_col++) {
			printf("Bbox for %d %d: top = %d, bottom = %d, left = %d, right = %d\n", c_row, c_col, row_loc[c_row], row_loc[c_row+1], col_loc[c_col], col_loc[c_col+1]);
		}

	}
*/

	for (c_row = 1; c_row < row_loc[0] - 1; c_row++) {
		for (c_col = 1; c_col < col_loc[0]; c_col++) {
			struct recognition *recog = recognize_char(image, max, row_loc[c_row], row_loc[c_row+1], col_loc[c_col], col_loc[c_col+1]);
			if (recog->guess[0].codepoint == 'X') {
				query_char(image, max, row_loc[c_row], row_loc[c_row+1], col_loc[c_col], col_loc[c_col+1]);
			}
			free(recog);
		}
		putchar('\n');

	}

	return 0;
}

/* Throw down row boundaries all in one go.
 */
int *rowdown(gray **image, gray max, int direction, int num_rows, int num_cols, float distance)
{
	float *counts;
	int *splits;
	float best_cost = INFINITY;
	float best_dist = 0, best_offset = 0;
	float current_dist;
	int dimension = (direction == HORZ ? num_rows : num_cols);
	int j;
	float k;

	counts = countup(image, max, direction, num_rows, num_cols);

	for (current_dist = distance - 2;
	     current_dist < distance + 2;
	     current_dist += 0.01) {
		float current_offset;
		for (current_offset = 0;
		     current_offset < distance + 2;
		     current_offset++) {
			float current_cost = 0;
			for (k = current_offset; k < dimension; k += current_dist) {
				current_cost += counts[(int)floor(k)];
			}
			if (current_cost < best_cost) {
				best_cost = current_cost;
				best_dist = current_dist;
				best_offset = current_offset;
			}

		}
	}

	splits = malloc(sizeof(int) * ceil(dimension / best_dist + 1));

	j = 1;
	for (k = best_offset; k < dimension; k += best_dist) {
		splits[j++] = (int)floor(k);
	}
	splits[0] = j-1;

	free(counts);
	return splits;

}

/* Count the values in each row in the image.
 * If given HORZ, find a "vertical" array as we are finding baselines.
 */
float *countup(gray **image, gray max, int direction, int num_rows, int num_cols)
{
	float *totals;
	int i;
	int dimension = (direction == HORZ ? num_rows : num_cols);

	totals = malloc(sizeof(float) * (dimension));

	for (i = 0; i < dimension; i++) {
		totals[i] = row_count(image, max, direction, num_rows, num_cols, i);
	}

	return totals;
}

/* `offset' is the row or column number to add up.
 * If given HORZ, count across the image; VERT count up and down.
 */
float row_count(gray **image, gray max, int direction, int num_rows, int num_cols, int offset)
{
	int x = 0;
	float sum = 0;

	if ((offset < 0) || (offset > (direction == HORZ ?
				       num_rows :
				       num_cols)))
		// Can't count pixels that aren't there!
		return INFINITY;

	for ( x = 0; x < (direction == HORZ ?
			  num_cols :
			  num_rows) ; x++ ) {
		float pixel = (float)(max - (direction == HORZ ?
					     image[offset][x] :
					     image[x][offset]))
			       / (float)max;
		sum += pixel * pixel;
	}

	return sum;
}


/* Display function for console training
 */
int print_char(gray **image, gray max, int top, int bottom, int left, int right)
{
	int row, col;

	int overshoot_y = 3;
	int overshoot_x = 10;

	char line[300];
	int i;
	char figures[] = {' ', '.', '*', '%', 'A', '#', '@'};

	for (row = top - overshoot_y; row < bottom + overshoot_y; row++) {
		i = 0;
		if ((row == top) || (row == bottom)) {
			for (col = left - overshoot_x ; col < right + overshoot_x; col++) {
				if ((col == left) || (col == right)) {
					line[i++] = '+';
				}
				line[i++] = '-';
				line[i++] = '-';
			}
			line[i++] = '\0';
			puts(line);
			i = 0;
		}
		for (col = left - overshoot_x; col < right + overshoot_x; col++) {
			if ((col == left) || (col == right)) {
				line[i++] = '|';
			}
			line[i++] = figures[(sizeof(figures)-1)-image[row][col]/(max/sizeof(figures))];
			line[i++] = figures[(sizeof(figures)-1)-image[row][col]/(max/sizeof(figures))];
		}
		line[i++] = '\0';
		puts(line);
	}
	return 0;
}

/* Ask the user to identify a character.
 */
int query_char(gray **image, gray max, int top, int bottom, int left, int right)
{
	char input[10];

	print_char(image, max, top, bottom, left, right);
	puts("What is this?");
	fgets(&input, 9, stdin);
	train_char(image, max, top, bottom, left, right, (int)input[0]);
}

/* Attempt to identify a character mechanically.
 */
struct recognition *recognize_char(gray **image, gray max, int top, int bottom, int left, int right)
{
	int row, col;
	gray value = 0;
	struct recognition *found = malloc(sizeof(struct recognition));

	for (row = top; row < bottom; row++) {
		for (col = left; col < right ; col++) {
			value += image[row][col];
		}
	}

	value /= (bottom-top) * (right-left);
	if (value > 192)
		found->guess[0].codepoint = ' ';
	else
		found->guess[0].codepoint = 'X';
}

/* Compare the small image `reference' to a region of the large image `test'.
 * Returns a estimate of similarity.
 */
float compare_char(gray **reference, int refwide, int refhigh, gray refmax, gray **test, gray testmax, int top, int bottom, int left, int right)
{
	int row, col;
	int error;
	int test_adjust;

	test_adjust = refmax / testmax;

	error = 0;
	for (row = 0; row < refhigh; row++)
		for (col = 0; col < refwide; col++) {
			error += abs(reference[row][col] - test[row + top][col + left] * test_adjust);
		}
	printf("error = %d\n", error);
	return (float)error;
}

/* Put a character into the set of trained pictures.
 */
int train_char(gray **image, gray max, int top, int bottom, int left, int right, int codepoint)
{
	int i = 0;
	int overshoot = 1;
	int t_row, t_col;
	struct trained_picture *chosen_one;

	if (training_data == NULL) {
		training_data = malloc(128 * sizeof(struct trained_char));
		training_size = 128;
		memset(training_data, 0, training_size * sizeof(struct trained_char));
	}

	while (1) {
		if (training_data[i].codepoint == UNDEFINED_CODEPOINT) {
			/* This is a character that is new to us */
			training_data[i].codepoint = codepoint;
			training_data[i].pictures = malloc(sizeof(struct trained_picture));
			chosen_one = training_data[i].pictures;
			chosen_one->data = pgm_allocarray(right-left + 2 * overshoot, bottom-top + 2 * overshoot);
			chosen_one->example_count = 0;
			chosen_one->next = NULL;

			for (t_row = 0; t_row < bottom-top + 2*overshoot; t_row++)
				for (t_col = 0; t_col < right-left + 2*overshoot; t_col++)
					chosen_one->data[t_row][t_col] = 0;

			break;
		} else if (training_data[i].codepoint == codepoint) {
			/* Hmm, seen one of these already */
			chosen_one = training_data[i].pictures;
			break;
		} else if (i > training_size) {
			/* Training table is full!
			 * XXX REALLOC
			 */
			return;
		} else {
			i++;
			continue;
		}
	}

	t_row = 0;
	t_col = 0;
	for (int row = top - overshoot;
	     row < bottom + overshoot;
	     row++) {
		t_col = 0;
		for (int col = left - overshoot;
		     col < right + overshoot;
		     col++) {
			chosen_one->data[t_row][t_col] += image[row][col];
			t_col++;
		}
		t_row++;
	}
	chosen_one->example_count++;
	return;
}

