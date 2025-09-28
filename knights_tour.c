#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 100
#define M N
#define HEAP_MAX_SIZE 8

#define UPDATE_ADJACENT_DEGREES(pos, delta)                                                        \
    do {                                                                                           \
        const int base_i = (pos).i;                                                                \
        const int base_j = (pos).j;                                                                \
        if (base_i + 2 < N && base_j + 1 < M && !visited[(base_i + 2) * M + base_j + 1])           \
            board[(base_i + 2) * M + base_j + 1] += (delta);                                       \
        if (base_i + 2 < N && base_j - 1 >= 0 && !visited[(base_i + 2) * M + base_j - 1])          \
            board[(base_i + 2) * M + base_j - 1] += (delta);                                       \
        if (base_i - 2 >= 0 && base_j + 1 < M && !visited[(base_i - 2) * M + base_j + 1])          \
            board[(base_i - 2) * M + base_j + 1] += (delta);                                       \
        if (base_i - 2 >= 0 && base_j - 1 >= 0 && !visited[(base_i - 2) * M + base_j - 1])         \
            board[(base_i - 2) * M + base_j - 1] += (delta);                                       \
        if (base_i + 1 < N && base_j + 2 < M && !visited[(base_i + 1) * M + base_j + 2])           \
            board[(base_i + 1) * M + base_j + 2] += (delta);                                       \
        if (base_i + 1 < N && base_j - 2 >= 0 && !visited[(base_i + 1) * M + base_j - 2])          \
            board[(base_i + 1) * M + base_j - 2] += (delta);                                       \
        if (base_i - 1 >= 0 && base_j + 2 < M && !visited[(base_i - 1) * M + base_j + 2])          \
            board[(base_i - 1) * M + base_j + 2] += (delta);                                       \
        if (base_i - 1 >= 0 && base_j - 2 >= 0 && !visited[(base_i - 1) * M + base_j - 2])         \
            board[(base_i - 1) * M + base_j - 2] += (delta);                                       \
    } while (0)

typedef struct {
    int i;
    int j;
} Coord;

typedef struct {
    unsigned degree : 4;
    unsigned visited : 1;
} Square;

typedef struct {
    Coord pos;
    int move_index;
} StackFrame;

typedef struct {
    StackFrame* frames;
    int top;
} Stack;

void initialize_stack(Stack* stack) {
    stack->frames = (StackFrame*)malloc(N * M * sizeof(StackFrame));
    stack->top = -1;
}

void free_stack(Stack* stack) {
    free(stack->frames);
}

void stack_push(Stack* stack, Coord pos, int move_index) {
    if (stack->top < N * M - 1) {
        stack->frames[++stack->top] = (StackFrame){pos, move_index};
    }
}

StackFrame* stack_top(Stack* stack) {
    return stack->top >= 0 ? &stack->frames[stack->top] : NULL;
}

void stack_pop(Stack* stack) {
    if (stack->top >= 0) {
        stack->top--;
    }
}

int stack_is_empty(Stack* stack) {
    return stack->top < 0;
}

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
    Coord* coords;
    int length;
} Path;

void initialize_path(Path* path) {
    path->coords = (Coord*)malloc(N * M * sizeof(Coord));
    path->length = 0;
}

void free_path(Path* path) {
    free(path->coords);
}

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
    { 2,  1},
    { 2, -1},
    {-2,  1},
    {-2, -1},
    { 1,  2},
    { 1, -2},
    {-1,  2},
    {-1, -2}
};

inline int is_valid_position(Coord c) {
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

inline int coord_flatten(Coord c) {
    return c.i * M + c.j;
}

int* initialize_board() {
    int* board = (int*)malloc(N * M * sizeof(int));
    Coord c;
    for (c.i = 0; c.i < N; ++c.i) {
        for (c.j = 0; c.j < M; ++c.j) {
            board[coord_flatten(c)] = calculate_degree(c);
        }
    }
    return board;
}

int* initialize_visited() {
    int* visited = (int*)malloc(N * M * sizeof(int));
    for (int i = 0; i < N * M; ++i) {
        visited[i] = 0;
    }
    return visited;
}

Square* initialize_squares() {
    Square* squares = (Square*)malloc(N * M * sizeof(Square));
    Coord c;
    for (c.i = 0; c.i < N; ++c.i) {
        for (c.j = 0; c.j < M; ++c.j) {
            squares[coord_flatten(c)].degree = calculate_degree(c);
            squares[coord_flatten(c)].visited = 0;
        }
    }
    return squares;
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

int solve_iterative(int* board, int* visited, Path* path, Coord start) {
    Stack stack;
    initialize_stack(&stack);

    path_push(path, start);
    visited[coord_flatten(start)] = 1;

    UPDATE_ADJACENT_DEGREES(start, -1);

    stack_push(&stack, start, -1);

    while (!stack_is_empty(&stack) && path->length < N * M) {
        StackFrame* frame = stack_top(&stack);
        const Coord current = frame->pos;

        MinHeap heap;
        init_heap(&heap);

        for (int k = 0; k < 8; ++k) {
            const Coord next = coord_add(current, KNIGHT_MOVES[k]);
            if (is_valid_position(next)) {
                const int idx = coord_flatten(next);
                if (!visited[idx]) {
                    const Move move = {.degree = board[idx], .move_idx = k, .coord = next};
                    heap_insert(&heap, move);
                }
            }
        }

        int found_move = 0;

        while (!heap_is_empty(&heap) && !found_move) {
            const Move move = heap_extract_min(&heap);
            const Coord next = move.coord;
            const int next_idx = coord_flatten(next);

            if (!visited[next_idx]) {
                path_push(path, next);
                stack_push(&stack, next, move.move_idx);
                visited[next_idx] = 1;
                UPDATE_ADJACENT_DEGREES(next, -1);
                found_move = 1;
            }
        }

        if (!found_move) {
            if (path->length > 1) {
                path_pop(path);
                visited[coord_flatten(current)] = 0;
                UPDATE_ADJACENT_DEGREES(current, 1);
            }
            stack_pop(&stack);
        }
    }

    free_stack(&stack);
    return path->length == N * M;
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

    printf("Starting knight's tour from (%d, %d) on a %dx%d board.\n", start.i, start.j, N, M);

    int* board = initialize_board();     // Need to free
    int* visited = initialize_visited(); // Need to free

    Path path;
    initialize_path(&path); // Need to free

    Coord c = start;

    clock_t t_start = clock();

    solve_iterative(board, visited, &path, c);

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

    free(board);
    free(visited);
    free_path(&path);

    return 0;
}
