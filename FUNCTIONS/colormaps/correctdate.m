%% From ExampleYearsTracking
function [datelist] = correctdate(datesin)
datelist = zeros(1,length(datesin));
for i = 1:length(datesin)
    date2=datestr(datesin(i));
    daynf=str2double(date2(1:2));
    month=date2(4:6);
    [monthint, ~] = monthconversion(month);
    year=str2double(date2(8:11));
    date3 = datetime([year monthint daynf]);
    datelist(i) = day(date3,'dayofyear');
    %center around the new year
    if (datelist(i) > 300)
        if (mod(i,4) == 0)
            leap = 366;
        else
            leap = 365;
        end
        datelist(i) = datelist(i)-leap;
    end
end
end