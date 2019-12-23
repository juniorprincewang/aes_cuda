all:
	nvcc aes.cu -o aes --cudart=shared

clean:
	rm -rf aes