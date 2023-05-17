ffmpeg -framerate 10 -start_number 501 -i density_proj_%05d.png -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2, fps=25,format=yuv420p" density_proj_destroy.mp4
