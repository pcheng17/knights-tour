#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 100
#define M N
#define HEAP_MAX_SIZE 8

typedef struct {
    int i;
    int j;
} Coord;

typedef struct {
    int degree;
    int move_idx;
    Coord coord;
} Move;

typedef struct {
    Move moves[HEAP_MAX_SIZE];
    int size;
} MinHeap;

void init_heap(MinHeap* heap) {
    heap->size = 0;
}

void swap_moves(Move* a, Move* b) {
    Move temp = *a;
    *a = *b;
    *b = temp;
}

void heapify_up(MinHeap* heap, int index) {
    while (index > 0) {
        const int parent = (index - 1) / 2;
        if (heap->moves[parent].degree > heap->moves[index].degree) {
            swap_moves(&heap->moves[index], &heap->moves[parent]);
            index = parent;
        } else {
            break;
        }
    }
}

void heapify_down(MinHeap* heap, int index) {
    while (1) {
        int smallest = index;
        const int left = 2 * index + 1;
        const int right = 2 * index + 2;

        if (left < heap->size && heap->moves[left].degree < heap->moves[smallest].degree) {
            smallest = left;
        }
        if (right < heap->size && heap->moves[right].degree < heap->moves[smallest].degree) {
            smallest = right;
        }
        if (smallest != index) {
            swap_moves(&heap->moves[index], &heap->moves[smallest]);
            index = smallest;
        } else {
            break;
        }
    }
}

void heap_insert(MinHeap* heap, Move move) {
    if (heap->size == HEAP_MAX_SIZE) return;
    heap->moves[heap->size] = move;
    heapify_up(heap, heap->size++);
}

Move heap_extract_min(MinHeap* heap) {
    Move minMove = heap->moves[0];
    heap->moves[0] = heap->moves[--heap->size];
    heapify_down(heap, 0);
    return minMove;
}

int heap_is_empty(MinHeap* heap) {
    return heap->size == 0;
}

Coord coord_add(Coord a, Coord b) {
    return (Coord){a.i + b.i, a.j + b.j};
}

typedef struct {
    Coord coords[N * M];
    int length;
} Path;

void path_push(Path* path, Coord pos) {
    if (path->length < N * M) {
        path->coords[path->length++] = pos;
    }
}

void path_pop(Path* path) {
    if (path->length > 0) {
        --path->length;
    }
}

Coord KNIGHT_MOVES[8] = {
    {2, 1}, {2, -1}, {-2, 1}, {-2, -1},
    {1, 2}, {1, -2}, {-1, 2}, {-1, -2}
};

int is_valid_position(Coord c) {
    return (c.i >= 0 && c.i < N && c.j >= 0 && c.j < M);
}

int calculate_degree(Coord c) {
    int d = 0;
    Coord next;
    for (int k = 0; k < 8; ++k) {
        next = coord_add(c, KNIGHT_MOVES[k]);
        if (is_valid_position(next)) {
            ++d;
        }
    }
    return d;
}

void fill_board_with_degrees(int board[N][M]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            board[i][j] = calculate_degree((Coord){i, j});
        }
    }
}

void print_board(int board[N][M]) {
    printf("Board of degrees (%dx%d):\n", N, M);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            printf("%d ", board[i][j]);
        }
        printf("\n");
    }
}

void initialize_visited(int visited[N][M]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            visited[i][j] = 0;
        }
    }
}

int get_max_degree(int board[N][M]) {
    int maxDegree = -1;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            if (board[i][j] > maxDegree) {
                maxDegree = board[i][j];
            }
        }
    }
    return maxDegree;
}

void write_tour_to_file(const char* filename, Path* path) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file for writing\n");
        return;
    }
    fprintf(file, "%d,%d\n", N, N); // First line is the size of the board.
    for (int i = 0; i < path->length; ++i) {
        fprintf(file, "%d,%d\n", path->coords[i].i, path->coords[i].j);
    }
    fclose(file);
}

void solve_recursive(int board[N][M], int visited[N][M], Path* path, Coord c, int max_degree) {
    path_push(path, c);
    visited[c.i][c.j] = 1;

    if (path->length == N * M) {
        return; // Tour complete
    }

    Coord next, adj;
    MinHeap heap;
    init_heap(&heap);

    for (int k = 0; k < 8; ++k) {
        next = coord_add(c, KNIGHT_MOVES[k]);
        if (is_valid_position(next) && !visited[next.i][next.j]) {
            Move move = {
                .degree = board[next.i][next.j],
                .move_idx = k,
                .coord = next
            };
            heap_insert(&heap, move);
        }
    }

    while (!heap_is_empty(&heap) && path->length < N * M) {
        Move move = heap_extract_min(&heap);
        next = move.coord;

        for (int k = 0; k < 8; ++k) {
            adj = coord_add(next, KNIGHT_MOVES[k]);
            if (is_valid_position(adj) && !visited[adj.i][adj.j]) {
                --board[adj.i][adj.j];
            }
        }

        solve_recursive(board, visited, path, next, max_degree);

        if (path->length == N * M) {
            break; // Tour complete
        }

        // Backtrack
        for (int k = 0; k < 8; ++k) {
            adj = coord_add(next, KNIGHT_MOVES[k]);
            if (is_valid_position(adj) && !visited[adj.i][adj.j]) {
                ++board[adj.i][adj.j];
            }
        }

        path_pop(path);
        visited[next.i][next.j] = 0;
    }

    // Final backtrack for my current position
    if (path->length < N * M) {
        path_pop(path);
        visited[c.i][c.j] = 0;
        for (int k = 0; k < 8; ++k) {
            next = coord_add(c, KNIGHT_MOVES[k]);
            if (is_valid_position(next) && !visited[next.i][next.j]) {
                ++board[next.i][next.j];
            }
        }
    }
}

int main(int argc, char* argv[]) {

    Coord start = {0, 0};

    if (argc == 3) {
        start.i = atoi(argv[1]);
        start.j = atoi(argv[2]);
        if (!is_valid_position(start)) {
            fprintf(stderr, "Error: Invalid starting position (%d, %d)\n", start.i, start.j);
            return 1;
        }
    }

    int board[N][M];
    int visited[N][M];

    Path path = {.length = 0};

    initialize_visited(visited);
    fill_board_with_degrees(board);
    const int max_deg = get_max_degree(board);

    Coord c = start;

    clock_t t_start = clock();

    solve_recursive(board, visited, &path, c, max_deg);

    clock_t t_end = clock();
    double time_spent = (double)(t_end - t_start) / CLOCKS_PER_SEC;

    if (path.length == N * M) {
        printf("Success! The knight's tour is complete.\n");
    } else {
        printf("The knight's tour ended prematurely after %d moves.\n", path.length);
        printf("The knight could not visit all %d squares.\n", N * M);
    }

    printf("Time taken: %f seconds\n", time_spent);

    write_tour_to_file("knights_tour.txt", &path);
    printf("Knight's tour path written to knights_tour.txt\n");

    return 0;
}
