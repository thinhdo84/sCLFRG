function F=normalized_form(S)
A=S.A;b=S.b;
F=[];
for i =1:size(A,1)
    F =  [F;A(i,:)/b(i)];
end

end