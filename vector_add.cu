#include <cuda_runtime.h>
#include <stdio.h>

__global__ void vectorAdd(int *a, int *b, int *c, int n) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i < n)
    c[i] = a[i] + b[i];
}

int main() {
  int n = 1000000;
  size_t size = n * sizeof(int);
  int *h_a, *h_b, *h_c;
  int *d_a, *d_b, *d_c;

  // Allocate host memory
  h_a = (int *)malloc(size);
  h_b = (int *)malloc(size);
  h_c = (int *)malloc(size);

  // Initialize host arrays
  for (int i = 0; i < n; i++) {
    h_a[i] = i;
    h_b[i] = i * 2;
  }

  // Allocate device memory
  cudaMalloc((void **)&d_a, size);
  cudaMalloc((void **)&d_b, size);
  cudaMalloc((void **)&d_c, size);

  // Copy input data from host to device memory
  cudaMemcpy(d_a, h_a, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, h_b, size, cudaMemcpyHostToDevice);

  // Launch kernel
  int threadsPerBlock = 256;
  int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
  vectorAdd<<<blocksPerGrid, threadsPerBlock>>>(d_a, d_b, d_c, n);

  // Copy result from device to host
  cudaMemcpy(h_c, d_c, size, cudaMemcpyDeviceToHost);

  // Verify result
  for (int i = 0; i < n; i++) {
    if (h_c[i] != h_a[i] + h_b[i]) {
      fprintf(stderr, "Result verification failed at element %d!\n", i);
      exit(1);
    }
  }

  printf("Test PASSED\n");

  // Free device memory
  cudaFree(d_a);
  cudaFree(d_b);
  cudaFree(d_c);

  // Free host memory
  free(h_a);
  free(h_b);
  free(h_c);

  return 0;
}