%Function to read data from the .txt file
%Opens the document and sets the 'r' as a reading doc
doc = fopen(filename,'r');

disp('|-------------------------------------|'); 
disp(['| Reading data from "', filename, '.txt" |']);
disp('|-------------------------------------|');
aux = false;
while feof(doc)~=1
    
    %Clears the text variables, since it changes during the document
    text = []; linsplit = [];
    
    %Reads the line of the document
    line = fgetl(doc);
    
    %Verifies if the line has texts
    if ~isempty(line)
        %Verifies if the line has separations
        if ~contains(line,':')
            text = line;  %If the line does not have the separation, then the text is the line
        else
            linsplit = strread(line,'%s','delimiter',':'); %If there is separation, the text is the first word
            text = char(linsplit(1));
        end 
    end
    
    if  ~isempty(line)
        if ~strcmp(text(1), '%')
            if strcmp(text, 'Body') %Identifies the beginnin of the Body header
                nbody = str2double(linsplit(2));
                aux = true;
            elseif aux == true
                info = strread(line,'%s','delimiter',':');
                if strcmp(info(1),'Name')
                    Bodies.B(nbody).(char(info(1))) = char(info(2));
                else
                    Bodies.B(nbody).(char(info(1))) = str2num(char(info(2)));
                end
            end
        end
    end
    
    if  ~isempty(line)
        if ~strcmp(text(1), '%')
            if strcmp(text, 'System') %Identifies the beginnin of the Body header
                
            elseif aux == false
                info = strread(line,'%s','delimiter',':');
                Bodies.System.(char(info(1))) = str2num(char(info(2)));
            end
        end
    end
    
end
clear aux;
disp(['Number of bodies: ', num2str(size(Bodies.B,2))]);
disp('Number of flexible modes: ');
nflex = [];
for i = 1:size(Bodies.B,2)
    disp(['Body ',num2str(i),' ( ',Bodies.B(i).Name,' ) : ',num2str(sum(Bodies.B(i).Flexible_modes))]);
    Bodies.B(i).nflex = sum(Bodies.B(i).Flexible_modes);
    nflex = [nflex Bodies.B(i).nflex];
end

for i = 1:size(Bodies.B,2)
   aux = 6*i+sum(nflex(1:i))-(6+nflex(i))+1;
   Bodies.B(i).pos_v = [aux:aux+2];
   Bodies.B(i).pos_omega = [aux+3:aux+5];
   if nflex(i) > 0
       Bodies.B(i).pos_qfp = [aux+6:aux+6+nflex(i)-1];
   else
       Bodies.B(i).pos_qfp = [aux+5:aux+5];
   end
   aux = [0,cumsum(nflex)];
   el = length(nflex)*6+sum(nflex)+3*(i-1)+aux(i)+1;
   Bodies.B(i).pos_R = [el:el+2];
   if nflex(i) > 0
       Bodies.B(i).pos_qf = [el+3:el+3+nflex(i)-1];
   else
       Bodies.B(i).pos_qf = [el+2:el+2];
   end
   el = 9*length(nflex)+sum(nflex)*2+1+(3+Bodies.System.quat)*(i-1);
   Bodies.B(i).pos_Theta = [el:el+2+Bodies.System.quat];
end

%To set the correct data type for some inputs
Bodies.System.Lock = logical(Bodies.System.Lock);
Bodies.System.Embedded = logical(Bodies.System.Embedded);