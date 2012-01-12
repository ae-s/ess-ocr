/* Split scans of listing files into character cells.
 * Duncan Smith
 * 2009, 2011
 *
 * gcc ess-ocr.c -ggdb3 -Wall -lnetpbm -lm -std=c99
 */

//#define _POSIX_SOURCE

#include <sys/unistd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

#include <pgm.h>

#define NO_ROW -1

// Looking for a division between lines
#define HORZ 0

// Looking for a division between characters
#define VERT 1

// when moving things around to search for a fit, move by up to this
// distance.
#define JITTER 4

int *rowdown	(gray **image, gray max, int direction, int num_rows, int num_cols, float distance);

float *countup	(gray **image, gray max, int direction, int num_rows, int num_cols);
float row_count	(gray **image, gray max, int direction, int num_rows, int num_cols, int offset);

int print_char	(gray **image, gray max, int top, int bottom, int left, int right);
int query_char	(gray **image, gray max, int top, int bottom, int left, int right);
struct recognition *recognize_char	(gray **image, gray max,
					 int top, int bottom, int left, int right);
void train_char	(gray **image, gray max,
		 int top, int bottom, int left, int right,
		 int codepoint, int overshoot, int dirty);

int find_square (gray **image, gray max,
		 int xmax, int ymax,
		 int *ul_x, int *ul_y,
		 int *ur_x, int *ur_y,
		 int *ll_x, int *ll_y,
		 int *lr_x, int *lr_y);
void crop_to_rect	(gray ***image, int *cols, int *rows, gray max);
void write_debug	(gray **image, gray max, int cols, int rows, int col[], int row[]);
int training_save	(void);
int train_load_char	(char *filename);

float compare_char(gray **reference,
		   gray refmax,
		   int refhigh, int refwide,
		   gray **test, gray testmax,
		   int top, int bottom, int left, int right);

float weight_border(gray **img, gray max, int top, int bottom, int left, int right);

#define MAX_RECOGNITION 16
struct recognition {
	struct recognized_char {
		int codepoint;
		float certainty;
	} guess[MAX_RECOGNITION];
};

struct trained_picture {
	int dirty;	// 0 if this picture also exists on disk.
	int rows;
	int cols;
	int max;
	gray **data;
	struct trained_picture *next;	// linked list
};

#define UNDEFINED_CODEPOINT 0
#define MAX_CODEPOINT MAXINT
struct trained_char {
	int codepoint;
	struct trained_picture *pictures;	// head of linked list
};

struct trained_char *training_data;
int training_size = 0;

int main(int argc, char *argv[])
{
	gray **image;
	FILE *fp;
	FILE *outfile;
	gray max;
	int prows, pcols;
	int *row_loc;
	int *col_loc;

	char *out_page;
	int out_page_ptr;

	int c_row, c_col;

	pgm_init(&argc, argv);

	fp = fopen(argv[1], "r");
	image = pgm_readpgm(fp, &pcols, &prows, &max);
	fclose(fp);

	training_load();

	puts("Cropping");
	crop_to_rect(&image, &pcols, &prows, max);
	printf("- is now %d cols by %d rows\n", pcols, prows);

#define DPI 300
	puts("Finding columns");
	col_loc = rowdown(image, max, VERT, prows, pcols, .06375 * DPI);
	puts("Finding rows");
	row_loc = rowdown(image, max, HORZ, prows, pcols, .10750 * DPI);

	write_debug(image, max, pcols, prows, col_loc, row_loc);

	printf("%d rows and %d columns\n", row_loc[0], col_loc[0]);

/*	for (c_row = 1; c_row < row_loc[0] - 1; c_row++) {
		for (c_col = 1; c_col < col_loc[0]; c_col++) {
			printf("Bbox for %d %d: top = %d, bottom = %d, left = %d, right = %d\n", c_row, c_col, row_loc[c_row], row_loc[c_row+1], col_loc[c_col], col_loc[c_col+1]);
		}

	}
*/

	out_page = malloc(1024 * 32 * sizeof(char));

	for (c_row = 1; c_row < row_loc[0] - 1; c_row++) {
		for (c_col = 1; c_col < col_loc[0]; c_col++) {
			char chr;

			int x, y, xend, yend;

			/*
			 * how far to offset the character cell from
			 * the given one in order to ensure no blobs
			 * sitting across the cell borders.
			 */
			int xoff, yoff;
			float best_cost = 0;
			int best_xoff, best_yoff;

			char bar_x[] = "|--------o--------|";
			char bar_y[] = "|----o----|";

			y    = row_loc[c_row];
			yend = row_loc[c_row + 1];
			x    = col_loc[c_col];
			xend = col_loc[c_col + 1];

			/* find a good place to drop the character cell border */
			for (xoff = -JITTER*2; xoff <= JITTER*2; xoff++)
				for (yoff = -JITTER; yoff <= JITTER; yoff++) {
					float cost;
					cost = weight_border(image, max, y + yoff, yend + yoff, x + xoff, xend + xoff);
					if (cost > best_cost) {
						best_xoff = xoff;
						best_yoff = yoff;
						best_cost = cost;
					}
				}

			bar_x[9 + best_xoff] = '@';
			bar_y[5 + best_yoff] = '@';

			printf("Jittering border by %s in x and %s in y\n", bar_x, bar_y);

			x    += best_xoff;
			xend += best_xoff;
			y    += best_yoff;
			yend += best_yoff;

			struct recognition *recog = recognize_char(image, max, y, yend, x, xend);
			if (recog != NULL) {
				chr = recog->guess[0].codepoint;
				if (chr == '\0') {
					chr = query_char(image, max, y, yend, x, xend);
				}
				free(recog);
			} else {
				chr = query_char(image, max, y, yend, x, xend);
			}
			out_page[out_page_ptr++] = chr;
			out_page[out_page_ptr + 1] = '\0';
			printf("\e[H\e[2J");

			//printf("\e[H");
			puts(out_page);
			printf("ptr %d\n", out_page_ptr);
		}
		out_page[out_page_ptr++] = '\n';
	}

	outfile = fopen(argv[2], "w");
	fwrite(out_page, sizeof(char), out_page_ptr, outfile);
	fclose(outfile);

	return 0;
}

float weight_border(gray **img, gray max, int top, int bottom, int left, int right)
{
	int off;
	float total = 0;
	float divisor = 0;

	for (off = top; off < bottom; off++) {
		total += img[off][left] + img[off][right];
		divisor += max * 2;
	}

	for (off = left + 1; off < right - 1; off++) {
		total += img[top][off] + img[bottom][off];
		divisor += max * 2;
	}

	return total / divisor;
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

	printf("Splits for direction %d ... \n", direction);
	printf("  suggested at intervals of %f\n", distance);
	printf("  are from %f at intervals of %f, with goodness of %f.\n", best_offset, best_dist, best_cost);

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
 *
 * returns (0) for full white row, (1 * row size) for full black row.
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
			line[i++] = figures[(sizeof(figures))-image[row][col]/(max/sizeof(figures))];
			line[i++] = figures[(sizeof(figures))-image[row][col]/(max/sizeof(figures))];
		}
		line[i++] = '\0';
		puts(line);
	}
	return 0;
}

int print_char_plain(gray **image, gray max, int rows, int cols)
{
	char figures[] = {' ', '.', '*', '%', 'A', '#', '@'};

	printf("This image is %d rows by %d cols up to %d\n", rows, cols, max);

	char line[cols + 2];
	for (int y = 0; y < rows; y++) {
		int x;
		for (x = 0; x < cols; x++) {
			line[x] = figures[(sizeof(figures))-image[y][x]/(max/sizeof(figures))];
			printf("%04x ", image[y][x]);
		}
		line[++x] = '\0';
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
	train_char(image, max, top, bottom, left, right, (int)input[0], 1, 1);
	training_save();
	return (int) input[0];
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
	if (value > 240) {
		// very white, probably a space
		found->guess[0].codepoint = ' ';
	} else {
		// not very white, probably a character
		struct trained_char *cur_test;
		float cur_best_val = INFINITY;
		struct trained_char *best_test;
		float best_val = INFINITY;

		int test_count = 0;

		cur_test = training_data;
		if (cur_test == NULL) return NULL;

		while (training_data[test_count].codepoint != UNDEFINED_CODEPOINT) {
			struct trained_picture *pic;
			cur_best_val = INFINITY;
			pic = cur_test->pictures;
			while (pic != NULL) {
				float bval;
				int xo, yo;
				for (xo = -JITTER; xo < JITTER; xo++)
					for (yo = -JITTER; yo < JITTER; yo++) {
						bval = compare_char(pic->data, pic->max,
								    pic->rows, pic->cols,
								    image, max,
								    top + yo, bottom + yo, left + xo, right + xo);
						if (bval < cur_best_val)
							cur_best_val = bval;
					}

				pic = pic->next;
			}

			if (cur_best_val < best_val) {
				best_val = cur_best_val;
				best_test = cur_test;
			}

			cur_test++;
			test_count++;
		}

		found->guess[0].codepoint = best_test->codepoint;
		found->guess[0].certainty = best_val;

		if (best_val > 0.1)
			found->guess[0].codepoint = '\0';
	}

	return found;
}

/* Compare the small image `reference' to a region of the large image `test'.
 * Returns a estimate of similarity.
 */
float compare_char(gray **reference,
		   gray refmax,
		   int refhigh, int refwide,
		   gray **test, gray testmax,
		   int top, int bottom, int left, int right)
{
	int row, col;
	float error;
	float test_adjust;

	test_adjust = refmax / testmax;
	error = 0;

	for (row = 0; row < refhigh; row++) {
		for (col = 0; col < refwide; col++) {
			gray pix_ref, pix_test;
			pix_ref = reference[row][col];
			pix_test = test[row + top][col + left];
			error += fabs((float)(pow(((float)pix_ref - (float)pix_test) * test_adjust, 2)));
		}
	}

	error = sqrt(error / (row * col)) / testmax;

//	printf("%d,%d / %d,%d - ", top, left, bottom, right);
//	printf("compares to %f\n", error);
	return error;
}

/* Put a character into the set of trained pictures.
 */
void train_char(gray **image, gray max, int top, int bottom, int left, int right, int codepoint, int overshoot, int dirty)
{
	int i = 0;
	int t_row, t_col;
	struct trained_picture *chosen_one;

	if (training_data == NULL) {
		training_size = 128;
		training_data = malloc(training_size * sizeof(struct trained_char));
		memset(training_data, 0, training_size * sizeof(struct trained_char));
	}

	chosen_one = malloc(sizeof(struct trained_picture));
	chosen_one->dirty = dirty;
	chosen_one->rows = bottom-top + 2 * overshoot;
	chosen_one->cols = right-left + 2 * overshoot;
	chosen_one->data = pgm_allocarray(chosen_one->cols, chosen_one->rows);
	chosen_one->max = max;
	chosen_one->next = NULL;

	for (t_row = 0; t_row < bottom-top + 2*overshoot; t_row++)
		for (t_col = 0; t_col < right-left + 2*overshoot; t_col++) {
			(chosen_one->data)[t_row][t_col] = image[top + t_row][left + t_col];
		}

	while (1) {
		if (training_data[i].codepoint == UNDEFINED_CODEPOINT) {
			/* This is a character that is new to us */
			puts("new character");
			training_data[i].codepoint = codepoint;
			training_data[i].pictures = chosen_one;

			break;
		} else if (training_data[i].codepoint == codepoint) {
			/* Hmm, seen one of these already */
			printf("adding new at %d\n", i);
			chosen_one->next = training_data[i].pictures;
			training_data[i].pictures = chosen_one;
			break;
		} else if (i > training_size) {
			puts("training table full!");
			/* Training table is full!
			 * XXX REALLOC
			 */
			return;
		} else {
			i++;
			continue;
		}
	}

/*
	t_row = 0;
	t_col = 0;
	for (int row = top - overshoot;
	     row < bottom + overshoot;
	     row++) {
		if (row > chosen_one->rows) break;
		t_col = 0;
		for (int col = left - overshoot;
		     col < right + overshoot;
		     col++) {
			if (col > chosen_one->cols) break;
			chosen_one->data[t_row][t_col] += image[row][col];
			t_col++;
		}
		t_row++;
	}
*/
	return;
}

#define SMALL_RATIO (0.2)

int find_square(gray **image, gray max_pix,
		int xmax, int ymax,
		int *ul_x, int *ul_y,
		int *ur_x, int *ur_y,
		int *ll_x, int *ll_y,
		int *lr_x, int *lr_y)
{
	int top = 0;
	int left = 0;
	int right = 0;
	int bot = 0;

	float max = INFINITY, peak = 0;
	int peak_i;

	printf("ymax %d, xmax %d\n", ymax, xmax);

	// 50 = skip black border
	for (int i = 50 ; i < (ymax * SMALL_RATIO) ; i++) {
		float sum = row_count(image, max_pix, HORZ, ymax, xmax, i);
		if (sum < max) {
			// typical white line
			max = sum;
		}
		if (sum > peak) {
			// found a border, probably
			peak = sum;
			peak_i = i;
		}

		// test for past-the-border
		if ((peak / max > (SMALL_RATIO)) && (fabs(sum / peak) < (1 - SMALL_RATIO)))
			top = peak_i;
	}

	max = INFINITY; peak = 0;
	// 50 = skip the black border around the page
	for (int i = ymax - 50 ; i > (ymax * (1 - SMALL_RATIO)) ; i--) {
		float sum = row_count(image, max_pix, HORZ, ymax, xmax, i);
		if (sum < max) {
			// found a typical white line
			max = sum;
		}
		if (sum > peak) {
			// found a border, probably
			peak = sum;
			peak_i = i;
		}

		// test for past-the-border
		if ((peak / max > SMALL_RATIO) && (fabs(sum / peak) < (1 - SMALL_RATIO)))
			bot = peak_i;
	}

	// 50 = skip black border
	max = INFINITY; peak = 0;
	for (int i = 50 ; i < (xmax * SMALL_RATIO) ; i++) {
		float sum = row_count(image, max_pix, VERT, ymax, xmax, i);
		if (sum < max) {
			// typical white line
			max = sum;
		}
		if (sum > peak) {
			// found a border, probably
			peak = sum;
			peak_i = i;
		}

		// test for past-the-border
		if ((peak / max > SMALL_RATIO) && (fabs(sum / peak) < (1 - SMALL_RATIO)))
			left = peak_i;
	}

	max = INFINITY; peak = 0;
	for (int i = xmax - 50 ; i > (xmax * (1 - SMALL_RATIO)) ; i--) {
		float sum = row_count(image, max_pix, VERT, ymax, xmax, i);
		if (sum < max) {
			// typical white line
			max = sum;
		}
		if (sum > peak) {
			// found a border, probably
			peak = sum;
			peak_i = i;
		}

		// test for past-the-border
		if ((peak / max > SMALL_RATIO) && (fabs(sum / peak) < (1 - SMALL_RATIO)))
			right = peak_i;
	}

	*ul_y = *ur_y = top;
	*ll_y = *lr_y = bot;
	*ul_x = *ll_x = left;
	*ur_x = *lr_x = right;

	return 0;
}

void crop_to_rect(gray ***image, int *cols, int *rows, gray max)
{
	int ul_x, ul_y, ur_x, ur_y, ll_x, ll_y, lr_x, lr_y;
	gray **img_new;
	int new_rows, new_cols;

	puts("squaring");
	find_square(*image, max, *cols, *rows,
		    &ul_x, &ul_y,
		    &ur_x, &ur_y,
		    &ll_x, &ll_y,
		    &lr_x, &lr_y);

	printf("from (%d, %d) to (%d, %d)\n", ul_x, ul_y, lr_x, lr_y);

	new_cols = (ur_x - ul_x);
	new_rows = (ll_y - ul_y);
	printf("cropping to %d cols by %d rows\n", new_cols, new_rows);
	img_new = pgm_allocarray(new_cols, new_rows);

	for (int x = ul_x; x < ur_x; x++) {
		for (int y = ul_y; y < ll_y; y++) {
			img_new[y - ul_y][x - ul_x] = (*image)[y][x];
		}
	}

	*rows = new_rows;
	*cols = new_cols;
	pgm_freearray(*image, *rows);
	*image = img_new;
}

/* create an image with dividing lines scribed on it
 */
void write_debug(gray **image, gray max, int cols, int rows, int col[], int row[])
{
	FILE *fp;
	gray *cur_row;
	gray *row_blacks;
	int row_i;

	row_blacks = malloc((cols) * sizeof(gray));
	cur_row = malloc((cols) * sizeof(gray));
	memset(row_blacks, 0, (cols) * sizeof(gray));

	fp = fopen("crop.pgm", "w");
	pgm_writepgminit(fp, cols, rows, max, 0);
	row_i = 1;

	for (int nrow = 0 ; nrow < rows; nrow++) {
		if (row[row_i] == nrow) {
			row_i++;
			pgm_writepgmrow(fp, row_blacks, cols, max, 0);
		} else {
			memcpy(cur_row, image[nrow], cols * sizeof(gray));
			for (int col_i = 1; col_i < col[0]; col_i++)
				cur_row[col[col_i]] = 0;
			pgm_writepgmrow(fp, cur_row, cols, max, 0);
		}
	}

	fclose(fp);
	return;
}

int training_save(void)
{
	struct trained_char *cur_point;

	cur_point = training_data;

	while (cur_point->codepoint != UNDEFINED_CODEPOINT) {
		struct trained_picture *cur_pic;
		cur_pic = cur_point->pictures;
		while (cur_pic != NULL) {
			if (cur_pic->dirty == 1) {
				char name[20];
				FILE *fp;

				sprintf(name, "./training/char_%04x_XXXXXX", cur_point->codepoint);
				mkstemp(name);

				fp = fopen(name, "w");
				pgm_writepgm(fp, cur_pic->data, cur_pic->cols, cur_pic->rows, cur_pic->max, 0);
				fclose(fp);

				cur_pic->dirty = 0;
			}
			cur_pic = cur_pic->next;
		}

		cur_point++;
	}

	return 0;
}

int training_load(void)
{
	DIR *entry;
	struct dirent *derp;
	int count = 0;
	struct stat st_buf;

	entry = opendir("./training");

	chdir("./training");

	derp = readdir(entry);
	while (derp != NULL) {
		stat(derp->d_name, &st_buf);
		if (S_ISREG(st_buf.st_mode))
			count += train_load_char(derp->d_name);

		derp = readdir(entry);
	}
	closedir(entry);

	chdir("..");

	return count;
}

int train_load_char(char *filename)
{
	gray **img;
	int rows, cols;
	gray max;
	FILE *fp;

	int chr;

	sscanf(filename, "char_%x", &chr);

	printf("reading %s into %02x(%c)\n", filename, chr, chr);

	fp = fopen(filename, "r");
	img = pgm_readpgm(fp, &cols, &rows, &max);
	fclose(fp);

	print_char_plain(img, max, rows, cols);

	train_char(img, max, 0, rows, 0, cols, chr, 0, 0);

	pgm_freearray(img, rows);
	return 1;
}
