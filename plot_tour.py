import argparse
from PIL import Image, ImageDraw

def draw_tour(positions, n, m):
    size = 80
    width, height = n * size, m * size
    background_color = "white"
    image = Image.new("RGB", (width+2, height+2), background_color)
    draw = ImageDraw.Draw(image)

    for i in range(n + 1):
        draw.line([(i * size, 0), (i * size, height)], fill= 'gray', width = 1)
    for i in range(m + 1):
        draw.line([(0, i * size), (width, i * size)], fill= 'gray', width = 1)

    path_coordinates = [(a * size + size/2, b* size+ size/2) for a, b in positions]
    line_width = 8
    L = len(path_coordinates)
    for i in range(L-1):
        line_color = "rgb(" + str(int(i/L * 255)) + ",0," + str(int(255 * (L-i)/L)) + ")"
        draw.line([path_coordinates[i], path_coordinates[i+1]], fill=line_color, width=line_width)

    return image

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("--file", type=str, default="knights_tour.txt", help="Path to the file containing the knight's tour.")
    arg_parser.add_argument("--output", type=str, help="Output image file name.")
    args = arg_parser.parse_args()

    with open("knights_tour.txt") as f:
        lines = f.readlines()
        n, m = map(int, lines[0].split(","))
        positions = [tuple(map(int, line.split(","))) for line in lines[1:]]

    if len(positions) != n * m:
        raise ValueError("The path length does not match the board size.")
    visited = set(positions)
    if len(visited) != n * m:
        raise ValueError("The path does not visit every cell exactly once.")
    for (x, y) in positions:
        if not (0 <= x < n and 0 <= y < m):
            raise ValueError(f"Position {(x, y)} is out of bounds.")

    print("All checks passed!")

    if args.output:
        start_pos = positions[0]
        idx = start_pos[0] * m + start_pos[1]
        image = draw_tour(positions, n, m)
        image.save(args.output)
        print(f"Image saved as {args.output}")
