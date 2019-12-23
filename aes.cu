#include "aes.h"

const char *file_path = "plaintext.txt";

//generate round keys from initial key
void expand_key(uint8_t *key, uint8_t *rkey){
  uint32_t i,j,k;
  uint8_t tempa[4];
  uint32_t nround = 10;

  //first round key is just the key
  for(i = 0; i < 4; ++i){
    rkey[4*i + 0] = key[4*i + 0];
    rkey[4*i + 1] = key[4*i + 1];
    rkey[4*i + 2] = key[4*i + 2];
    rkey[4*i + 3] = key[4*i + 3];
  }

  for(i = 4; i < 4*(nround + 1); ++i){
    for(j = 0; j < 4; ++j){
      tempa[j] = rkey[(i-1)*4 + j];
    }

    if(i % 4 == 0){
      //rotate 4 bytes in word
      k = tempa[0];
      tempa[0] = tempa[1];
      tempa[1] = tempa[2];
      tempa[2] = tempa[3];
      tempa[3] = k;


      tempa[0] = sbox[tempa[0]];
      tempa[1] = sbox[tempa[1]];
      tempa[2] = sbox[tempa[2]];
      tempa[3] = sbox[tempa[3]];
  
      tempa[0] = tempa[0] ^ rcon[i/4];

    }

    rkey[4*i + 0] = rkey[4*(i-4) + 0] ^ tempa[0];
    rkey[4*i + 1] = rkey[4*(i-4) + 1] ^ tempa[1];
    rkey[4*i + 2] = rkey[4*(i-4) + 2] ^ tempa[2];
    rkey[4*i + 3] = rkey[4*(i-4) + 3] ^ tempa[3];

  } 

}

//XOR round key with block(1 block per thread)
__device__ void add_round_key(uint8_t *block, uint8_t *key, uint32_t offset){
  //word size traversal
  uint32_t *b = (uint32_t *)block;
  uint32_t *k = (uint32_t *)key;
  for(int i = 0; i < 4; ++i){
    b[offset/4 + i] = b[offset/4 + i] ^ k[i];
  }  
}

//substitute block int sbox (1 block per thread)
__device__ void sub_bytes(uint8_t *block, uint32_t offset){
  for(int i = 0; i < 16; ++i){
    block[offset + i] = dsbox[block[offset + i]];
  }
}

//substitute block int sbox (1 block per thread)
__device__ void inv_sub_bytes(uint8_t *block, uint32_t offset){
  for(int i = 0; i < 16; ++i){
    block[offset + i] = disbox[block[offset + i]];
  }
}


//mix columns by taking linear combinations in the field (1 block per thread)
__device__ void mix_columns(uint8_t *block, uint32_t offset){
  for(int i = 0; i < 4; ++i){ //iterate over columns
    uint8_t a[4];
    uint8_t b[4]; 
    uint8_t h;
  
    for(int j = 0; j < 4; ++j){
      a[j] = block[offset + 4*i + j];
      h = (uint8_t)((int8_t)a[j] >> 7);
      b[j] = a[j] << 1;
      b[j] ^= 0x1b & h;
    } 

    block[offset + 4*i + 0] = b[0] ^ a[3] ^ a[2] ^ b[1] ^ a[1];
    block[offset + 4*i + 1] = b[1] ^ a[0] ^ a[3] ^ b[2] ^ a[2];
    block[offset + 4*i + 2] = b[2] ^ a[1] ^ a[0] ^ b[3] ^ a[3];
    block[offset + 4*i + 3] = b[3] ^ a[2] ^ a[1] ^ b[0] ^ a[0]; 

  }
}

//mix columns by taking linear combinations in the field (1 block per thread)
/**
f(x) = 11 * x^3 + 13 * x^2 + 9 * x +  14
*/
__device__ void inv_mix_columns(uint8_t *block, uint32_t offset){

	for(int i = 0; i < 4; ++i){ //iterate over columns
		uint16_t t;
		uint8_t a[4];
		uint8_t a2[4];
		uint8_t a4[4];
		uint8_t a8[4];
		a[0] = block[offset + 4*i];
		a[1] = block[offset + 4*i + 1];
		a[2] = block[offset + 4*i + 2];
		a[3] = block[offset + 4*i + 3];
		a2[0] = gmul2(a[0]);
		a2[1] = gmul2(a[1]);
		a2[2] = gmul2(a[2]);
		a2[3] = gmul2(a[3]);
		a4[0] = gmul4(a[0]);
		a4[1] = gmul4(a[1]);
		a4[2] = gmul4(a[2]);
		a4[3] = gmul4(a[3]);
		a8[0] = gmul8(a[0]);
		a8[1] = gmul8(a[1]);
		a8[2] = gmul8(a[2]);
		a8[3] = gmul8(a[3]);

		block[offset + 4*i + 0] = a8[0] ^ a4[0] ^ a2[0] ^ a8[1] ^ a2[1] ^ a[1] ^ a8[2] ^ a4[2] ^ a[2] ^ a8[3] ^ a[3];
		block[offset + 4*i + 1] = a8[0] ^ a[0] ^ a8[1] ^ a4[1] ^ a2[1] ^ a8[2] ^ a2[2] ^ a[2] ^ a8[3] ^ a4[3] ^ a[3];
		block[offset + 4*i + 2] = a8[0] ^ a4[0] ^ a[0] ^ a8[1] ^ a[1] ^ a8[2] ^ a4[2] ^ a2[2] ^ a8[3] ^ a2[3] ^ a[3];
		block[offset + 4*i + 3] = a8[0] ^ a2[0] ^ a[0] ^ a8[1] ^ a4[1] ^ a[1] ^ a8[2] ^ a[2] ^ a8[3] ^ a4[3] ^ a2[3];
	}
}


//shift rows left by 0,1,2,3 bytes respectively (1 block per thread)
__device__ void shift_rows(uint8_t *sblock, uint32_t offset){
  uint8_t tmp;

  uint8_t *block = sblock + offset; 

  //row 0 remains unshifted

  //shift row 1 left by 1
  tmp = block[1];
  block[1] = block[5];
  block[5] = block[9];
  block[9] = block[13];
  block[13] = tmp;

  //shift row 2 letf by 2
  tmp = block[2];
  block[2] = block[10];
  block[10] = tmp;

  tmp = block[6];
  block[6] = block[14];
  block[14] = tmp;

  //shift row 3 left by 3
  tmp = block[3];
  block[3] = block[15];
  block[15] = block[11];
  block[11] = block[7];
  block[7] = tmp;
}


//shift rows right by 0,1,2,3 bytes respectively (1 block per thread)
__device__ void inv_shift_rows(uint8_t *sblock, uint32_t offset){
  uint8_t tmp;
  uint8_t *block = sblock + offset; 

  //row 0 remains unshifted

  //shift row 1 right by 1
  tmp = block[13];
  block[13] = block[9];
  block[9] 	= block[5];
  block[5] = block[1];
  block[1] = tmp;

  //shift row 2 right by 2
  tmp = block[10];
  block[10] = block[2];
  block[2] = tmp;

  tmp = block[14];
  block[14] = block[6];
  block[6] = tmp;

  //shift row 3 right by 3
  tmp = block[3];
  block[3] = block[7];
  block[7] = block[11];
  block[11] = block[15];
  block[15] = tmp;
}

//aes 128 encryption with expanded key supplied
//implemented as basic byte algorithm (naive)
//operates on one block per thread
__device__ void encrypt(uint8_t *block, uint8_t *rkey, uint32_t offset){
	add_round_key(block, rkey, offset);
	for(int i = 1; i < 10; ++i){
		sub_bytes(block, offset);
		shift_rows(block, offset);
		mix_columns(block, offset);
		add_round_key(block, rkey + 16*i, offset);
	}
	sub_bytes(block, offset);
	shift_rows(block, offset);
	add_round_key(block, rkey + 160, offset);
}

__global__ void encrypt_one_block(uint8_t *block, uint8_t *rkey, uint32_t numblock)
{
	int bindex = blockIdx.x * blockDim.x + threadIdx.x;
	int offset = bindex * 16;
	if(bindex >= numblock) return;
	// printf("bindex = %d\n", bindex);
	encrypt(block, rkey, offset);
}

/***********************************************************************/
/* decrypt */
/***********************************************************************/

//aes 128 encryption with expanded key supplied
//implemented as basic byte algorithm (naive)
//operates on one block per thread
__device__ void decrypt(uint8_t *block, uint8_t *rkey, uint32_t offset){
	add_round_key(block, rkey+ 160, offset);
	for(int i = 9; i >=1 ; --i){
		inv_shift_rows(block, offset);
		inv_sub_bytes(block, offset);
		add_round_key(block, rkey + 16*i, offset);
		inv_mix_columns(block, offset);
	}
	inv_shift_rows(block, offset);
	inv_sub_bytes(block, offset);
	add_round_key(block, rkey, offset);
}

__global__ void decrypt_one_block(uint8_t *block, uint8_t *rkey, uint32_t numblock)
{
	int bindex = blockIdx.x * blockDim.x + threadIdx.x;
	int offset = bindex * 16;
	if(bindex >= numblock) return;
	// printf("bindex = %d\n", bindex);
	decrypt(block, rkey, offset);
}

void encrypt_cuda(uint8_t *data, uint8_t *out_data, uint8_t *key, uint32_t size)
{
	uint32_t numblock = size / 16;
	uint32_t num_bytes = size;
	uint8_t rkey[176];
	uint32_t *ddata;
	uint32_t *drkey;

	expand_key(key, rkey);
	cudaMalloc(&ddata, sizeof(uint8_t) * num_bytes);
	cudaMalloc(&drkey, sizeof(uint8_t) * 176);
	cudaMemcpy(ddata, (uint32_t *)data, sizeof(uint8_t) * num_bytes, cudaMemcpyHostToDevice);
	cudaMemcpy(drkey, (uint32_t *)rkey, sizeof(uint8_t) * 176, cudaMemcpyHostToDevice);


	encrypt_one_block<<<(numblock + 31)/32, 32>>>((uint8_t *)ddata, (uint8_t *)drkey, numblock);
	cudaThreadSynchronize();

	cudaMemcpy(out_data, ddata, sizeof(uint8_t) * num_bytes, cudaMemcpyDeviceToHost);
	//check for errors
	cudaError_t errCode = cudaPeekAtLastError();
	if(errCode != cudaSuccess){
	fprintf(stderr, "WARNING: A CUDA error occured: code=%d, %s\n", errCode, cudaGetErrorString(errCode));
	}
	cudaFree(ddata);
	cudaFree(drkey);
}

void decrypt_cuda(uint8_t *data, uint8_t *out_data, uint8_t *key, uint32_t size)
{
	uint32_t numblock = size / 16;
	uint32_t num_bytes = size;
	uint8_t rkey[176];
	uint32_t *ddata;
	uint32_t *drkey;

	expand_key(key, rkey);
	cudaMalloc(&ddata, sizeof(uint8_t) * num_bytes);
	cudaMalloc(&drkey, sizeof(uint8_t) * 176);
	cudaMemcpy(ddata, (uint32_t *)data, sizeof(uint8_t) * num_bytes, cudaMemcpyHostToDevice);
	cudaMemcpy(drkey, (uint32_t *)rkey, sizeof(uint8_t) * 176, cudaMemcpyHostToDevice);


	decrypt_one_block<<<(numblock + 31)/32, 32>>>((uint8_t *)ddata, (uint8_t *)drkey, numblock);
	cudaThreadSynchronize();

	cudaMemcpy(out_data, ddata, sizeof(uint8_t) * num_bytes, cudaMemcpyDeviceToHost);
	//check for errors
	cudaError_t errCode = cudaPeekAtLastError();
	if(errCode != cudaSuccess){
	fprintf(stderr, "WARNING: A CUDA error occured: code=%d, %s\n", errCode, cudaGetErrorString(errCode));
	}
	cudaFree(ddata);
	cudaFree(drkey);
}

void print_block_hex(const uint8_t *text)
{
	for(int i = 0; i < 16; ++i){
		printf("0x%x ", text[i]);
	}
	printf("\n");
}

void write_file(char *file)
{
	uint8_t plaintext[] = { 0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee ,0xff
		, 0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee ,0xff};
	FILE *fp;
	fp = fopen(file, "wb+");
	if(!fp) {
		fprintf(stderr, "FAILED: to open plaintext!\n");
		exit(-1);
	}
	fwrite(plaintext, 1, sizeof(plaintext), fp);
	fclose(fp);
}

//read data from file
uint8_t *file_buf(const char *file, long int *size){
	int fd = open(file, O_RDONLY);
	struct stat stats;
	if(fd < 0) {
		fprintf(stderr, "Error opening file\n");
		exit(1);
	}
  	if(fstat(fd, &stats) < 0) {
  		fprintf(stderr, "Error opening file\n");
  		exit(1);
  	}
  	uint8_t *mem = (uint8_t *)mmap(NULL, stats.st_size, PROT_READ, MAP_PRIVATE, fd, 0); 
  	if(mem == MAP_FAILED) {
  		fprintf(stderr, "mmap failed\n");
  		exit(1);
  	}
  	*size = stats.st_size;
  	return mem;
}


/*
Plaintext:  00112233445566778899aabbccddeeff
Cipher key: 000102030405060708090a0b0c0d0e0f
Ciphertext: 69c4e0d86a7b0430d8cdb78070b4c55a
*/
void aes_128_testcase_single()
{
	uint8_t plaintext[16] = { 0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee ,0xff};
	uint8_t ciphertext[16] = {0x69, 0xc4, 0xe0, 0xd8, 0x6a, 0x7b, 0x04, 0x30, 0xd8, 0xcd, 0xb7, 0x80, 0x70, 0xb4, 0xc5, 0x5a};
	uint8_t key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f};
	long int size = sizeof(plaintext);
	uint8_t *output_gpu = (uint8_t *) malloc(sizeof(uint8_t) * size);

	printf("Plaintext: \n");
	for(int i=0; i<size/16; i++) {
		print_block_hex(plaintext+16*i);
	}
	
	encrypt_cuda(plaintext, output_gpu, key, size);
	printf("Output: \n");
	for(int i=0; i<size/16; i++) {
		print_block_hex(output_gpu+16*i);
	}

	decrypt_cuda(ciphertext, output_gpu, key, size);
	printf("Output: \n");
	for(int i=0; i<size/16; i++) {
		print_block_hex(output_gpu+16*i);
	}
}

void aes_128_testcase_file()
{
	long int size;
	uint8_t key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f};
	uint8_t *plaintext = file_buf(file_path, &size);
	uint8_t *output_gpu = (uint8_t *) malloc(sizeof(uint8_t) * size);
	
	printf("Plaintext: \n");
	for(int i=0; i<size/16; i++) {
		print_block_hex(plaintext+16*i);
	}
	
	encrypt_cuda(plaintext, output_gpu, key, size);
	printf("Output: \n");
	for(int i=0; i<size/16; i++) {
		print_block_hex(output_gpu+16*i);
	}

	decrypt_cuda(output_gpu, output_gpu, key, size);
	printf("Output: \n");
	for(int i=0; i<size/16; i++) {
		print_block_hex(output_gpu+16*i);
	}
}

void aes_128_testcase_1M()
{
	long int size;
	char *file = "myfile1M";
	uint8_t key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f};
	uint8_t *plaintext = file_buf(file, &size);
	uint8_t *output_gpu = (uint8_t *) malloc(sizeof(uint8_t) * size);
	
	printf("Plaintext first 16 bytes: \n");
	print_block_hex(plaintext);
	
	encrypt_cuda(plaintext, output_gpu, key, size);

	decrypt_cuda(output_gpu, output_gpu, key, size);
	printf("Output first 16 bytes: \n");
	print_block_hex(output_gpu);
}

void aes_128_testcase_128M()
{
	long int size;
	char *file = "myfile128M";
	uint8_t key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f};
	uint8_t *plaintext = file_buf(file, &size);
	uint8_t *output_gpu = (uint8_t *) malloc(sizeof(uint8_t) * size);
	
	printf("Plaintext first 16 bytes: \n");
	print_block_hex(plaintext);
	
	encrypt_cuda(plaintext, output_gpu, key, size);

	decrypt_cuda(output_gpu, output_gpu, key, size);
	printf("Output first 16 bytes: \n");
	print_block_hex(output_gpu);
}


int main()
{
	// write_file(file_path);
	// return 0;
	/*
	uint8_t key[16] = {0x7E, 0x24, 0x06, 0x78, 0x17, 0xFA, 0xE0, 0xD7, 0x43, 0xD6, 0xCE, 0x1F, 0x32, 0x53, 0x91, 0x63};
	uint8_t rseed[16] = {0x00, 0x6C, 0xB6, 0xDB, 0xC0, 0x54, 0x3B, 0x59, 0xDA, 0x48, 0xD9, 0x0B, 0, 0, 0, 0};
	uint8_t plaintext[33] = "\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0A\x0B\x0C\x0D\x0E\x0F\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1A\x1B\x1C\x1D\x1E\x1F";
	long int size = 32;
	*/

	// aes_128_testcase_single();
	aes_128_testcase_128M();


	/*//Test correctness
	int sum = 0;
	for(uint32_t i = 0; i < size; ++i){
		sum += abs(ciphertext[i] - output_gpu[i]);
	}
	printf("Sum = %d  (should be zero if correct)\n", sum);*/
	return 0;
}