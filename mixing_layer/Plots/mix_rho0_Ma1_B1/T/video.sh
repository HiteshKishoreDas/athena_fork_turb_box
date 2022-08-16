ffmpeg -framerate 10 -i T_%05d.png -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2, fps=25, format=yuv420p" T_By.mp4
