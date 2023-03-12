#include <stdio.h>
#include <stdlib.h>

float** init_zero_vector(int m, int n){
  float* values = calloc(m*n, sizeof(float));
  float** vec = malloc(m*sizeof(float*));
  
  for (int i=0; i<m; ++i){
    vec[i] = values + i*n;
  }
  return vec;
}


void copy_vector(float** newV, float** previousV, int m, int n){
  size_t size = sizeof(float)*n;
  for (int i=0; i<m; ++i){
    for (int j=0; j<n; ++j){
      newV[i][j] = previousV[i][j];
    }
  }
}


void print_vector(float** vec, int m, int n){
  for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            printf("%f,", vec[i][j]);
        printf("\n");
    }
}


int write_result(char* fileName, float** vec, int m, int n){
  FILE *outFile;
  outFile=fopen(fileName,"w+");
  if (outFile==NULL){
        printf("[Error] Failed to open output file %s", fileName);
        return -1;
    }

  for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            fprintf(outFile, "%f,", vec[i][j]);
        fprintf(outFile, "\n");
    }
  fclose(outFile);

  return 0;
}