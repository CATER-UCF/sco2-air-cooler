import os
from PIL import Image


def main(folder):
    for file in os.listdir(folder):
        if file.endswith('.png') or file.endswith('.PNG'):
            file_path = os.path.join(folder, file)
            img = Image.open(file_path)
            img.save(file_path.replace('.png', '.tiff')
                     .replace('.PNG', '.tiff'))


if __name__ == '__main__':
    main('C:\\Users\\Rusty\\Downloads\\sCO2_Air_Cooler')
