%% Month Function to convert strings ("Jan") to other strings ("01") and ints (01)
function [monthint, monthintstr] = monthconversion(monthstr)
if monthstr == "Jan"
    monthintstr = "01";
    monthint = 01;
elseif monthstr == "Feb"
    monthintstr = "02";
    monthint = 02;
elseif monthstr == "Mar"
    monthintstr = "03";
    monthint = 03;
elseif monthstr == "Apr"
    monthintstr = "04";
    monthint = 04;
elseif monthstr == "May"
    monthintstr = "05";
    monthint = 05;
elseif monthstr == "Jun"
    monthintstr = "06";
    monthint = 06;
elseif monthstr == "Jul"
    monthintstr = "07";
    monthint = 07;
elseif monthstr == "Aug"
    monthintstr = "08";
    monthint = 08;
elseif monthstr == "Sep"
    monthintstr = "09";
    monthint = 09;
elseif monthstr == "Oct"
    monthintstr = "10";
    monthint = 10;
elseif monthstr == "Nov"
    monthintstr = "11";
    monthint = 11;
elseif monthstr == "Dec"
    monthintstr = "12";
    monthint = 12;
end
end