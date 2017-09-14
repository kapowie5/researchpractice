function thing = maketrap(x,a,b,c,d)
%a,b,c,d are the points of the trapezoid that will go from 0 to 1, with
%length d-a=len
    len=length(x);
    thing = zeros(1,len);
    for i = a:b
        thing(i) = 1/(b-a)*(i-a);
    end
    for i = b:c
        thing(i) = 1;
    end
    for i = c:d
        thing(i) = abs(1/(d-c)*i + 1/(c-d)*c - 1);
    end
    
%     https://www.mathworks.com/help/fuzzy/trapmf.html
%     thing = zeros(1,len);
%     for i = a:b
%         thing(i) = (x-a)/(b-a);
%     end
%     for i = b:c
%         thing(i) = 1;
%     end
%     for i = c:d
%         thing(i) = (d-x)/(d-c);
%     end

end
    