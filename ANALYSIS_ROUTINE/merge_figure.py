from PIL import Image
import argparse
import numpy as np

def create_canvas(index_list, folder_path, output_filename, row, column):
    index = 1
    total_width = 0
    total_height = 0
    for i in range(0, row):
        for j in range(0, column):
            filename = folder_path + "rho_x_0_0000" + str(index) + ".jpg"
            img = Image.open(filename)
            width, height = img.size
            index = index + 1
            if i == 0: total_width = total_width + width
        total_height = total_height + height


    merge_image = Image.new(mode="RGB", size=(total_width, total_height), color=(255,255,255))
    return merge_image


def main(index_list, folder_path, output_filename, row, column):

    merge_image = create_canvas(index_list, folder_path, output_filename, row, column)

    index = 1
    for i in range(0, row):
        for j in range(0, column):
            filename = folder_path + "rho_x_0_0000" + str(index) + ".jpg"
            img = Image.open(filename)
            width, height = img.size
            merge_image.paste(img, (width * j, height * i))
            index = index + 1

    merge_image.save(folder_path + "test.pdf")

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--row", help = "number of row in grid")
    parser.add_argument("-c", "--column", help = "number of column in grid")
    parser.add_argument("-i", "--input", help = "input image index")
    parser.add_argument("-f", "--folder", help= "root file folder path")
    parser.add_argument("-o", "--output", help= "output fileaname")
    args = parser.parse_args()


    if args.row != "": 
        row = int(args.row)

    if args.column != "":
        column = int(args.column)

    if args.input != "":
        index_list = args.input.split(",")

        if len(index_list) != (row * column):
            raise Exception("len of index list is %d. Not compatiable in %d x %d grid" %(len(index_list), row, column))

    if args.output != "":
        output_filename = args.output

    if args.folder != "":
        folder_path = args.folder

    

     
    main(index_list, folder_path, output_filename, row, column)
