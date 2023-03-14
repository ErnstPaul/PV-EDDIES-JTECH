try
echo on
mex inpoly.c
echo off
catch exception
    if(~isempty(exception))
        fprintf(['\n Error during compilation, be sure :\n'...
            'i)  You have C compiler installed (prefered compiler are MSVC/Intel/GCC)\n'...
            'ii) You did "mex -setup" in matlab prompt before running mexme_inpoly\n']);
    end    
end