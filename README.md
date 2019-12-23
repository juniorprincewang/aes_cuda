# aes_cuda
Simple AES-128 encryption/decryption implementation on CUDA.

To generate random file:  
```
head -c 1M </dev/urandom >myfile1M
```

Or if `head` doesn't understand the `M` suffix, then specify the size in bytes:  
```
head -c 1048576 </dev/urandom >myfile1M
```
