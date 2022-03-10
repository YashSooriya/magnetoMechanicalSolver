function dbfact=dbfact(i)

dbfact = 1;
if i>0
    for j=i:-1:1 %[i,1,-1]
        dbfact=dbfact*j;
    end
elseif i==0
    dbfact=1;
else
    error(message('error in factoral'));
end