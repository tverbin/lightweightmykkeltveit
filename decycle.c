
#include <inttypes.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#define BASE_BITS (2)
#define BASE (1l << BASE_BITS)
#define MAX_DIGIT (BASE - 1)
#define PRECISION (0.0001f)


enum mykkeltveit_color {
	MC_C = 0,
	MC_L,
	MC_R,
};

void print_num_in_base_4(int64_t num, int64_t digits)
{
	for (int64_t i = digits - 1; i >= 0; i--)
		printf("%ld", (num >> (i * BASE_BITS)) & MAX_DIGIT);
	printf("\n");
}

enum mykkeltveit_color get_mykkeltveit_color(int64_t kmer, int64_t k, float *x_values_arr)
{
	float sum = 0;
	int i;

	for (i = 0; i < k; i++) {
		sum += x_values_arr[k-i-1] * (float)(((MAX_DIGIT << (i*BASE_BITS)) & kmer) >> (i*BASE_BITS));
	}
	if (sum > -PRECISION && sum < PRECISION)
		return MC_C;
	else if (sum > 0)
		return MC_L;
	else
		return MC_R;
}

int64_t choose_vertex_to_remove(int64_t vertex, int64_t k, int64_t num_of_k_mers, float *x_values_arr)
{
	int64_t last_char;
	int64_t current_vertex = vertex;
	int64_t next_vertex;

	enum mykkeltveit_color current_color;
	enum mykkeltveit_color next_color;


	last_char = current_vertex >> (k * BASE_BITS - BASE_BITS);
	next_vertex = ((current_vertex & ((num_of_k_mers - 1) >> BASE_BITS)) << BASE_BITS) + last_char;

	current_color = get_mykkeltveit_color(vertex, k, x_values_arr);
	next_color = get_mykkeltveit_color(next_vertex, k, x_values_arr);

	/* TODO: replace to binary search implementation, possible due to the 
	 *       R-block L-block structure */
	while (next_color != MC_L || current_color == MC_L) {
		current_vertex = next_vertex;
		last_char = current_vertex >> (k * BASE_BITS - BASE_BITS);
		next_vertex = ((current_vertex & ((num_of_k_mers - 1) >> BASE_BITS)) << BASE_BITS) + last_char;

		if (next_vertex == vertex) { /* all vertices in cycle are color 0 */
			return vertex;
		}

		current_color = next_color;
		next_color = get_mykkeltveit_color(next_vertex, k, x_values_arr);
	}

	return next_vertex;
}

void print_decycle_set(int64_t k)
{
	int64_t num_of_k_mers = 1l << (k * BASE_BITS); // 4^k = 2^2^k
	float *x_values_arr; // the weight of each digit in PCR
	int64_t cycle = 0;
	int64_t j;
	int64_t i;
	int64_t remove_vertex;
	int64_t c_j_i;
	int64_t count = 0;

	x_values_arr = malloc(sizeof(float) * k);
	if (x_values_arr == NULL) {
		perror("malloc");
		exit(1);
	}

	for(i = 0; i < k; i++) {
		x_values_arr[i] = sinf(2 * (float)M_PI * i / k);
	}

	/* FKM algorithm implementation */
	i = k;
	do {
		if (k % i == 0) {
			count++;
			/* the vertex is symmetric and as a result the whole will be 
			 * colored with C */
			if (i != k)
				remove_vertex  = cycle;
			/* choose vertex based on mykkeltveit de-Bruijn coloring */
			else
				remove_vertex = choose_vertex_to_remove(cycle, k, 
						num_of_k_mers, x_values_arr);
			print_num_in_base_4(remove_vertex, k);
		}

		i = k;

		while ((((cycle >> (BASE_BITS * (i-1))) & MAX_DIGIT) == MAX_DIGIT) && (i != 0))
			i--;
		cycle += (1l << (BASE_BITS * (i-1)));

		for (j = i; j < k; j++) {
			c_j_i = (cycle & (MAX_DIGIT << (BASE_BITS*(j-i)))) << (BASE_BITS*i);
			cycle = (cycle & ((1l << (BASE_BITS*k)) - 1 - (MAX_DIGIT << (BASE_BITS*j)))) + c_j_i;
		}
		if (i == 0l) {
			break;
		}
	} while (i != 0);

	free(x_values_arr);

	fprintf(stderr, "total %ld %ld-mers from %ld\n", count, k, num_of_k_mers);
}

int main(int argc, char *argv[])
{
	int64_t k;
	if (argc < 2) {
		fprintf(stderr, "Usage: %s <k>\n", argv[0]);
		exit(1);
	}

	k = strtol(argv[1], NULL, 0);

	print_decycle_set(k);
	
	return 0;
}

