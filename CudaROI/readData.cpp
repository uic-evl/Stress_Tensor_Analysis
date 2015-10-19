#include <cstdlib>
#include <fstream>
#include <iostream>

int main(int argc, char** argv){
	
	FILE *fp = fopen(argv[1], "rb");
    if (!fp) {
        fprintf(stderr, "Error opening file '%s'\n", argv[1]);
        return 0;
    }

	size_t size = 32*32*32;
    void *data = malloc(size);
    size_t read = fread(data, 1, size, fp);
    fclose(fp);

	char* output = "out.txt";
	
	FILE *fout = fopen(output, "w+");
	size_t write = fwrite ( data, size, 1, fout );
	fclose (fout);

    return 0;
}