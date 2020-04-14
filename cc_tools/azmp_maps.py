# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 13:50:49 2019

@author: GIBBO
"""
import pandas as pd
import cc_tools as cct
import os

my_file = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_plot.xlsx')
years = ['2014', '2015', '2016', '2017', '2018']
seasons = ['spring','summer','fall']
variables = ['pH','Omega_A']
depth = ['surface', 'bottom']
for yy in years:
    for season in seasons:
        for d in depth:
            for var in variables:
                print yy, season, d, var
                cct.azmp_map(my_file, yy, season, d, var)



import numpy as np
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None

os.chdir('C:\Users\gibbo\Documents\carbonates\AZMP_OA\AZMP_OA data\AZMP_OA plots\AZMP_OA maps')
images_list = os.listdir('C:\Users\gibbo\Documents\carbonates\AZMP_OA\AZMP_OA data\AZMP_OA plots\AZMP_OA maps')
years = ['2014', '2015', '2016', '2017', '2018'] 
seasons = ['spring','summer','fall']
count = 0
for year in years:
    for season in seasons:
        background = Image.new('RGBA', (5438, 4500), (255, 255, 255, 255))
        imgs_s = [item for item in images_list if year in item and season in item and 'surface' in item]          
        imgs_b = [item for item in images_list if year in item and season in item and 'bottom' in item]          
        imgs_s = [Image.open(i) for i in imgs_s]
        imgs_b = [Image.open(i) for i in imgs_b]
        if  imgs_s == []:
            print ('no data that season')
            continue
        count += 1
        min_img_shape_s = sorted ([(np.sum(i.size), i.size) for i in imgs_s])[0][1]
        min_img_shape_b = sorted ([(np.sum(i.size), i.size) for i in imgs_b])[0][1]
        img_merge_s = np.hstack((np.asarray(i.resize(min_img_shape_s,Image.ANTIALIAS)) for i in imgs_s))
        img_merge_s = Image.fromarray(img_merge_s).save('AZMP_OA_s.png')
        img_merge_b = np.hstack((np.asarray(i.resize(min_img_shape_b,Image.ANTIALIAS)) for i in imgs_b))
        img_merge_b = Image.fromarray(img_merge_b).save('AZMP_OA_b.png')
        images = ['AZMP_OA_s.png', 'AZMP_OA_b.png']
        imgs = [Image.open(i) for i in images]
        min_img_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
        combined = np.vstack((np.asarray(i.resize(min_img_shape,Image.ANTIALIAS)) for i in imgs)) 
        combined = Image.fromarray(combined).save('AZMP_OA_' + year + '_' + season +'.png')
        ##########image final size w=5438, h=4150###############        
        image = Image.open('AZMP_OA_' + year + '_' + season +'.png').convert('RGB')
        font = ImageFont.truetype('arial.ttf', 100)
        draw = ImageDraw.Draw(image).text((250,100), 'a)', (255,255,255), font=font)
        draw = ImageDraw.Draw(image).text((2975,100), 'b)', (255,255,255), font=font)
        draw = ImageDraw.Draw(image).text((250,2230), 'c)', (255,255,255), font=font)
        draw = ImageDraw.Draw(image).text((2975,2230), 'd)', (255,255,255), font=font)
        #draw = ImageDraw.Draw(image).rectangle(((2400,200), (3200,380)), fill='black')
        #draw = ImageDraw.Draw(image).text((2550, 230), (season + ' ' + year), font=font)    
        cap_font = ImageFont.truetype('arial.ttf', 100, encoding='unic')
        #draw = ImageDraw.Draw(background).text((50, 4200), u'\u03a9', (0,0,0), font=cap_font)
        draw = ImageDraw.Draw(background).text((50, 4200), ('Figure ' + str(count) + ': Bathymetric map of the Atlantic Zone during '+ season + ' ' + year+', for ' + u'\u03a9'+ 'arg and pH at surface (a, b) and bottom (c, d).'), (0,0,0), font=cap_font)    
        background.paste(image, (0,0))
        image_final = background.save('C:\Users\gibbo\Documents\carbonates\AZMP_OA\AZMP_OA data\AZMP_OA plots\AZMP_OA_combinedmaps\AZMP_OA_' + year + '_' + season +'.png')     
        
        
os.chdir('C:\Users\gibbo\Documents\carbonates\AZMP_OA\AZMP_OA data\AZMP_OA plots\AZMP_OA_combinedmaps')        
images_list = os.listdir('C:\Users\gibbo\Documents\carbonates\AZMP_OA\AZMP_OA data\AZMP_OA plots\AZMP_OA_combinedmaps')
images_list.sort(key=os.path.getmtime)
figs = [item for item in images_list]          
figs = [Image.open(i).convert('RGB') for i in figs]
pdf_filename = ('C:\Users\gibbo\Documents\carbonates\AZMP_OA\AZMP_OA data\AZMP_OA plots\AZMP_OA_allfigs.pdf')
figs[0].save(pdf_filename, 'PDF', save_all=True, append_images=figs[1:])


