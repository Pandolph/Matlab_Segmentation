function new(S)
center=[62,71]; 
m=[55;67;68;39];
n=[72;74;34;68];
SegOutline = zeros(100);
newSegOutline = SegOutline;
for i = 1:length(m(:))
    if m(i)>=center(1)
        if n(i)>=center(2)
            newSegOutline(m(i):m(i)+2,n(i):n(i)+2)=1;
        else
            newSegOutline(m(i):m(i)+2,n(i)-2:n(i))=1;
        end
    else        
        if n(i)>=center(2)
            newSegOutline(m(i)-2:m(i),n(i):n(i)+2)=1;
        else
            newSegOutline(m(i)-2:m(i),n(i)-2:n(i))=1;
        end
    end
end
newSegOutline = newSegOutline-SegOutline;