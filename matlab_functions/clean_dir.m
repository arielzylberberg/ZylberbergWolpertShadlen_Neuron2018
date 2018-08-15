function clean_dir()
%borra los archivos con "~" del directorio actual
if isunix
    !rm *~*
else
    !del *~
end