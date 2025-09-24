#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 150

int KNIGHT_MOVES[8][2] = {
    {2, 1}, {2, -1}, {-2, 1}, {-2, -1},
    {1, 2}, {1, -2}, {-1, 2}, {-1, -2}
};

int is_valid_position(int i, int j) {
    return (i >= 0 && i < N && j >= 0 && j < N);
}

int calculate_degree(int i, int j) {
    int d = 0;
    for (int k = 0; k < 8; ++k) {
        int iNext = i + KNIGHT_MOVES[k][0];
        int jNext = j + KNIGHT_MOVES[k][1];

        if (is_valid_position(iNext, jNext)) {
            ++d;
        }
    }
    return d;
}

void fill_board_with_degrees(int board[N][N]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            board[i][j] = calculate_degree(i, j);
        }
    }
}

void print_board(int board[N][N]) {
    printf("Board of degrees (%dx%d):\n", N, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            printf("%d ", board[i][j]);
        }
        printf("\n");
    }
}

void initialize_visited(int visited[N][N]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            visited[i][j] = 0;
        }
    }
}

int get_max_degree(int board[N][N]) {
    int maxDegree = -1;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (board[i][j] > maxDegree) {
                maxDegree = board[i][j];
            }
        }
    }
    return maxDegree;
}

int main(int argc, char* argv[]) {

    int iStart = 0, jStart = 0; // Default starting position
    if (argc == 3) {
        iStart = atoi(argv[1]);
        jStart = atoi(argv[2]);
        if (!is_valid_position(iStart, jStart)) {
            fprintf(stderr, "Error: Invalid starting position (%d, %d)\n", iStart, jStart);
            return 1;
        }
    }

    int board[N][N];
    int visited[N][N];

    FILE* file = fopen("knights_tour.txt", "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file for writing\n");
        return 1;
    }

    fprintf(file, "%d,%d\n", N, N); // First line is the size of the board.

    initialize_visited(visited);
    fill_board_with_degrees(board);
    const int max_deg = get_max_degree(board);

    int i = iStart, j = jStart;
    visited[i][j] = 1;
    int moves = 1;

    clock_t start = clock();

    while (1) {
        fprintf(file, "%d,%d\n", i, j);
        int minDegreeMoveIdx = -1, minDegree = max_deg + 1;
        for (int k = 0; k < 8; ++k) {
            int iNext = i + KNIGHT_MOVES[k][0];
            int jNext = j + KNIGHT_MOVES[k][1];

            if (is_valid_position(iNext, jNext) && !visited[iNext][jNext]) {
                int degree = board[iNext][jNext];
                if (degree < minDegree) {
                    minDegreeMoveIdx = k;
                    minDegree = degree;
                }
            }
        }

        if (minDegreeMoveIdx == -1) {
            break; // No more moves possible
        }

        i += KNIGHT_MOVES[minDegreeMoveIdx][0];
        j += KNIGHT_MOVES[minDegreeMoveIdx][1];
        visited[i][j] = 1;
        ++moves;

        // Everything reachable from this new position needs to have its degree reduced by 1
        for (int k = 0; k < 8; ++k) {
            int iNext = i + KNIGHT_MOVES[k][0];
            int jNext = j + KNIGHT_MOVES[k][1];
            if (is_valid_position(iNext, jNext) && !visited[iNext][jNext]) {
                --board[iNext][jNext];
            }
        }
    }

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    if (moves == N * N) {
        printf("Success! The knight's tour is complete.\n");
    } else {
        printf("The knight's tour ended prematurely after %d moves.\n", moves);
        printf("The knight could not visit all %d squares.\n", N * N);
    }

    printf("Time taken: %f seconds\n", time_spent);

    fclose(file);
    printf("Knight's tour path written to knights_tour.txt\n");

    return 0;
}
