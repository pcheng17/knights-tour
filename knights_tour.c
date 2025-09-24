#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 6
#define M N

typedef struct {
    int i;
    int j;
} Coord;

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

    Coord next;
    int minDegreeMoveIdx = -1, minDegree = max_degree + 1;

    for (int k = 0; k < 8; ++k) {
        next = coord_add(c, KNIGHT_MOVES[k]);
        if (is_valid_position(next) && !visited[next.i][next.j]) {
            const int degree = board[next.i][next.j];
            if (degree < minDegree) {
                minDegreeMoveIdx = k;
                minDegree = degree;
            }
        }
    }

    if (minDegreeMoveIdx != -1) {
        next = coord_add(c, KNIGHT_MOVES[minDegreeMoveIdx]);

        Coord adj;
        for (int k = 0; k < 8; ++k) {
            adj = coord_add(next, KNIGHT_MOVES[k]);
            if (is_valid_position(adj) && !visited[adj.i][adj.j]) {
                --board[adj.i][adj.j];
            }
        }

        solve_recursive(board, visited, path, next, max_degree);
    }

    // Backtrack
    if (path->length < N * M) {
        path_pop(path);
        visited[next.i][next.j] = 0;

        Coord adj;
        for (int k = 0; k < 8; ++k) {
            adj = coord_add(c, KNIGHT_MOVES[k]);
            if (is_valid_position(adj) && !visited[adj.i][adj.j]) {
                ++board[adj.i][adj.j];
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
