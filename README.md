1)	Instalar o compilador: https://www.msys2.org/
2)	Abir MSYS2 MINGW64 e colar o endereço da pasta de destino: 
Ex: User@DESKTOP-U7KH33E MINGW64 //wsl.localhost/Ubuntu-24.04/home/gabrielsilverio/Doutorado/Disciplinas/USDM_model
3)	Comando para compilar em Fortran (rodar do terminal MSYS2 MINGW64):
gfortran USRMOD.f90 -o usermod64.dll -shared -fno-underscoring -static
