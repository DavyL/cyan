add_executable ( fourier1D fourier1D.c)
target_link_libraries ( fourier1D cyan )

add_test (NAME fourier1D COMMAND fourier1D)

