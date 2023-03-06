from PIL import Image
import argparse
import numpy as np


def create_canvas(input_filename, row, column):
    img = Image.open(input_filename)
    width, height = img.size
    total_width = column * width
    total_height = row * height
    merge_image = Image.new(mode="RGB", size=(total_width, total_height), color=(255, 255, 255))
    return merge_image


def main(folder_path, output_filename, row, column, var):
    input_filename = folder_path + var + "_x_0_00001" + ".jpg"
    merge_image = create_canvas(input_filename, row, column)

    index = 0
    for i in range(0, row):
        for j in range(0, column):
            filename = folder_path + "rho_x_0_0000" + str(index_list[index]) + ".jpg"
            img = Image.open(filename)
            width, height = img.size
            merge_image.paste(img, (width * j, height * i))
            index = index + 1

    merge_image.save(folder_path + output_filename)

    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--row", help="number of row in grid")
    parser.add_argument("-c", "--column", help="number of column in grid")
    parser.add_argument("-i", "--input", help="input image index")
    parser.add_argument("-f", "--folder", help="root file folder path")
    parser.add_argument("-o", "--output", help="output fileaname")
    parser.add_argument("-v", "--variables", help="variable to merge")
    args = parser.parse_args()

    if args.row != "":
        row = int(args.row)
    else:
        raise Exception("row must > 0")

    if args.column != "":
        column = int(args.column)
    else:
        raise Exception("column must > 0")

    if args.input != "":
        index_list = args.input.split(",")

        if len(index_list) != (row * column):
            raise Exception("len of index list is %d. Not compatiable in %d x %d grid" % (len(index_list), row, column))
    else:
        raise Exception("must be specific the index list")

    if args.output != "":
        output_filename = args.output
    else:
        raise Exception("must be specific the output filename")

    if args.folder != "":
        folder_path = args.folder
    else:
        raise Exception("must be specific the root folder")

    if args.variables:
        var = args.variables
    else:
        raise Exception("must be specific the variables")

    main(index_list, folder_path, output_filename, row, column, var)
