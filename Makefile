all:
	nvcc aes.cu -o aes

clean:
	rm -rf aes