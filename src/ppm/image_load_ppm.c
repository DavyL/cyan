#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <color/color.h>
#include <image/image.h>
#include <io/image_io.h>

#include "../config.h"

typedef struct {
	FILE *handle;
	char *line;
	size_t allocated_mem;	// taille allouée pour le tableau line
	ssize_t line_length;	// nbre de caractères contenus dans line ( incluant \n )
	int first_token;
} ppm_file_t;

ppm_file_t *ppm_file_new(char *filename);
void ppm_file_free(ppm_file_t * fichier);
ssize_t ppm_file_read_line(ppm_file_t * fichier);
char *ppm_file_get_next_token(ppm_file_t * fichier);
image_t *image_load_P3(ppm_file_t *, int *result);
image_t *image_load_P6(ppm_file_t *, int *result);


ppm_file_t *ppm_file_new(char *filename) {
	ppm_file_t *fichier = malloc(sizeof(ppm_file_t));
	fichier->handle = fopen(filename, "r");
	if (fichier->handle == (FILE *) NULL) {
		return NULL;
	}
	fichier->line = (char *) NULL;
	fichier->allocated_mem = (size_t) 0;
	fichier->line_length = (ssize_t) 0;
	fichier->first_token = 1;
	return fichier;
}

void ppm_file_free(ppm_file_t * fichier) {
	fclose(fichier->handle);
	if (fichier->line != (char *) NULL)
		free(fichier->line);
	free(fichier);
}

ssize_t ppm_file_read_line(ppm_file_t * fichier) {
	fichier->first_token = 1;
	fichier->line_length =
	    getline(&(fichier->line), &(fichier->allocated_mem),
		    fichier->handle);
	if (fichier->line_length > 0)
		if (fichier->line[0] == '#')
			return ppm_file_read_line(fichier);
	return (int) fichier->line_length;
}

char *ppm_file_get_next_token(ppm_file_t * fichier) {
	char delimiters[4] = { ' ', '\t', '\n', '\0' };
	char *chaine;
	char *result;
	if (fichier->line == NULL)
		ppm_file_read_line(fichier);
	if (fichier->first_token) {
		chaine = fichier->line;
		fichier->first_token = 0;
	} else
		chaine = NULL;
	result = strtok(chaine, delimiters);
	if (result == NULL) {
		ppm_file_read_line(fichier);
		result = ppm_file_get_next_token(fichier);
	}
	return result;
}


image_t *image_load_P3(ppm_file_t * inputfile, int *result) {

	image_t *img = NULL;
	int i, j;
	int width = atoi(ppm_file_get_next_token(inputfile));
	int height = atoi(ppm_file_get_next_token(inputfile));
	double maxval = (double) atoi(ppm_file_get_next_token(inputfile));

	fprintf(stderr, "P3 file:\n");
	fprintf(stderr, "\tsize: %d x %d \n", width, height);

	img = image_new(width, height, COLORTYPE_RGB);
	for (j = 0; j < img->rows; j++) {
		for (i = 0; i < img->cols; i++) {
			img->pixels[j * img->cols + i].coords[0] =
			    atoi(ppm_file_get_next_token(inputfile)) /
			    maxval;
			img->pixels[j * img->cols + i].coords[1] =
			    atoi(ppm_file_get_next_token(inputfile)) /
			    maxval;
			img->pixels[j * img->cols + i].coords[2] =
			    atoi(ppm_file_get_next_token(inputfile)) /
			    maxval;
		}
	}
	return img;
}


image_t *image_load_P6(ppm_file_t * inputfile, int *result) {

	int width = atoi(ppm_file_get_next_token(inputfile));
	int height = atoi(ppm_file_get_next_token(inputfile));
	int maxval = atoi(ppm_file_get_next_token(inputfile));

	fprintf(stderr, "P6 file:\n");
	fprintf(stderr, "\tsize: %d x %d \n", width, height);




	return NULL;		// FIXME 

}

image_t *image_load_ppm(char *filename, int *result) {
	image_t *img;
	ppm_file_t *inputfile;
	inputfile = ppm_file_new(filename);
	if (inputfile == NULL) {
		fprintf(stderr, "Cannot open file\n");
	}
	char *token;
	token = ppm_file_get_next_token(inputfile);

	if (strcmp(token, "P3") == 0) {
		img = image_load_P3(inputfile, result);
		*result = 1;
	} else {

		if (strcmp(token, "P6") == 0) {
			img = image_load_P6(inputfile, result);
			*result = 1;
		} else {
			fprintf(stderr, "Fichier non reconnu \n");
			*result = 0;
			img = NULL;
		}
	}
	ppm_file_free(inputfile);
	return img;
}