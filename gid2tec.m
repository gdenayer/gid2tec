function gid2tec(basename)

%% Name of carat or empire mesh and result files
nameFileMshInfo=strcat(basename,'.msh');
nameFileResInfo=strcat(basename,'.res');

if ~(isfile(nameFileMshInfo))
    fprintf('\t>> Error: The file %s does not exist.\n',nameFileMshInfo);
    return;
end
if ~(isfile(nameFileResInfo))
    fprintf('\t>> Error: The file %s does not exist.\n',nameFileResInfo);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract infos from the headers of the mesh
fid = fopen(nameFileMshInfo);
headersMsh=fgetl(fid);
fclose(fid);
headersMsh_info = strsplit(headersMsh);
if strcmp(headersMsh_info{1},'MESH')
    fprintf('\t>> The file %s is a MESH file...\n',nameFileMshInfo);
else
    fprintf('\t>> Error: The file %s is not a MESH file.\n',nameFileMshInfo);
    return;
end
headersMsh_info_tmp = strsplit(headersMsh,'dimension');
headersMsh_name=headersMsh_info_tmp{1};
headersMsh_info = strsplit(headersMsh_info_tmp{2});
headersMsh_info(1)= [];
headersMsh_dimension=str2num(headersMsh_info{1});
headersMsh_elementType=headersMsh_info{3};
headersMsh_Nnode=str2num(headersMsh_info{5});

% read the mesh and result files
fileMshInfo = fileread(nameFileMshInfo);
fileResInfo = fileread(nameFileResInfo);

%% Get the nodes and the elements of the meshes
%
% Get the nodes of the mesh
% the first stringSeparator is the header of the file
stringSeparator = headersMsh;
block = regexp(fileMshInfo,stringSeparator,'split');
% delete first index of 'block'
block(1) = [];
block = regexp(block{1},'Coordinates','split');
% delete first index of 'block'
block(1) = [];
% delete last index of 'block'
block(end) = [];
block = regexp(block{1},'End','split');
block(end) = [];
out = cell(size(block));
%
fprintf('\t>> Reading the coordinates of the nodes...');
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
msh.nodes = cell2mat(out);
msh.nodes = [msh.nodes(:,2:end) msh.nodes(:,1)]; % swap node ID to the last column
num_nodes = length(msh.nodes);
fprintf('DONE\n');
% tecplot seems to have some problem with the nodeID of FASTEST.
% The following reattributes an nodeID, which corresponds to the line in
% the msh file.
fprintf('\t>> Reformating the nodeIDs...');
maxNodeID=max(msh.nodes(:,4));
msh.nodeID2nodeLine=zeros(maxNodeID,1);
for i=1:num_nodes
    msh.nodeID2nodeLine(msh.nodes(i,4))=i;
end
fprintf('DONE\n');
%
% Get the elements of the mesh
fprintf('\t>> Reading the elements of the mesh...\n');
stringSeparator = 'Elements';
block = regexp(fileMshInfo,stringSeparator,'split');
block(1) = [];
block = regexp(block{1},'End','split');
block(end) = [];
out = cell(size(block));

if headersMsh_Nnode == 3
    % Triangle
    fprintf('\t>>   The MESH file contains triangles.\n');
    for k = 1:numel(block)
        out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1); % five values
        out{k} = horzcat(out{k}{:});
    end
elseif headersMsh_Nnode == 4
    % Quad
    fprintf('\t>>   The MESH file contains quad.\n');
    for k = 1:numel(block)
        out{k} = textscan(block{k},'%f %f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1); % five values
        out{k} = horzcat(out{k}{:});
    end
else
    fprintf('\t>> Error: ElementType not supported.\n');
    return;
end
elements = cell2mat(out);
% carat msh file introduces NaN in the elements array. This  is a trick to
% remove the lines with NaN
elements=elements(sum(isnan(elements),2)==0,:);
msh.elements = elements(:,2:end);
num_elements = length(msh.elements);
fprintf('\t>> DONE\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove header (it contains the word 'Result')
stringSeparator = 'GiD Post Results File 1.0';
block = regexp(fileResInfo,stringSeparator,'split');
% delete first index of 'block'
block(1) = [];
% create a cell array containing one result section per element
stringSeparator = 'Result';
block = regexp(block{1},stringSeparator,'split');
% delete first index of 'block'
block(1) = [];
nressection=size(block,2);
% look for the maximum number of step
resblock = regexp(block{end},'Values','split');
% extract info about the quantity in the current result section
% The delimiter has to be " for the case with variables including spaces
resblock_info=strsplit(resblock{1},'"');
resblock_info(1)=[];
outputQuantity_name=resblock_info{1};
% remove " from the name of quantity
outputQuantity_name=erase(outputQuantity_name,'"');
% replace whitespace with underscore
outputQuantity_name=regexprep(outputQuantity_name, ' ', '_');
resblock_info2=strsplit(resblock_info{4});
resblock_info2(1)=[];
ntotstep=str2num(resblock_info2{1})
% initialize for Tecplot file
nbquantityperstep=0;
outputQuantity_step_old=0;
% loop on all result sections
for ires=1:nressection
    fprintf('\t>> Reading results of the Result section number = %d/%d\n',ires,nressection);
    resblock = regexp(block{ires},'Values','split');
    % extract info about the quantity in the current result section
    % The delimiter has to be " for the case with variables including spaces
    resblock_info=strsplit(resblock{1},'"');
    resblock_info(1)=[];
    outputQuantity_name=resblock_info{1};
    % remove " from the name of quantity
    outputQuantity_name=erase(outputQuantity_name,'"');
    % replace whitespace with underscore
    outputQuantity_name=regexprep(outputQuantity_name, ' ', '_');
    resblock_info2=strsplit(resblock_info{4});
    resblock_info2(1)=[];
    outputQuantity_step=str2num(resblock_info2{1});
    outputQuantity_type=resblock_info2{2};
    outputQuantity_location=resblock_info2{3};
    if strcmpi(outputQuantity_type,'vector')
        %vector
        outputQuantity_name_x=strcat('"',outputQuantity_name,'_X','"');
        outputQuantity_name_y=strcat('"',outputQuantity_name,'_Y','"');
        outputQuantity_name_z=strcat('"',outputQuantity_name,'_Z','"');
        outputQuantity_varname=strcat(outputQuantity_name_x," ",outputQuantity_name_y," ",outputQuantity_name_z);
        if strcmp(outputQuantity_location,'OnNodes')
            outputQuantity_varlocation='NODAL NODAL NODAL';
        elseif strcmp(outputQuantity_location,'OnGaussPoints')
            outputQuantity_varlocation='CELLCENTERED CELLCENTERED CELLCENTERED';
        else
            fprintf('\t>> Error: Location not supported.\n');
            return;
        end
    else
        % scalar
        outputQuantity_varname=strcat('"',outputQuantity_name,'"');
        if strcmp(outputQuantity_location,'OnNodes')
            outputQuantity_varlocation='NODAL';
        elseif strcmp(outputQuantity_location,'OnGaussPoints')
            outputQuantity_varlocation='CELLCENTERED';
        else
            fprintf('\t>> Error: Location not supported.\n');
            return;
        end
    end
    %
    % write one tecplot file per step; A tecplot file may contain several
    % output quantities
    if outputQuantity_step > outputQuantity_step_old
        fprintf('\t>>   A new step starts.\n');
        if (nbquantityperstep>0)
            %% Write tecplot file
            fprintf('\t>>   Flush the previous data in a Tecplot file....\n');
            writeData2plt(basename,nameFileMshInfo,outputQuantity_step_old, ...
                outputQuantity,headersMsh_Nnode,num_nodes, ...
                num_elements,msh,nbquantityperstep,ntotstep)
        end
        nbquantityperstep=1;
    else
        nbquantityperstep=nbquantityperstep+1;
    end
    % store the infos for the current step
    outputQuantity.name(nbquantityperstep)=string(outputQuantity_name);
    outputQuantity.varname(nbquantityperstep)=string(outputQuantity_varname);
    outputQuantity.type(nbquantityperstep)=string(outputQuantity_type);
    outputQuantity.location(nbquantityperstep)=string(outputQuantity_location);
    outputQuantity.varlocation(nbquantityperstep)=string(outputQuantity_varlocation);
    %
    % remove the first element of the block: Typically the Result line
    resblock(1) = [];
    % remove the last elements of the block: Typically empty.
    resblock(2:end) = [];
    %resblock
    % Remove all after End
    resblock = regexp(resblock{1},'End','split');
    % Remove the last element of the block: Typically empty
    resblock(end) = [];
    out = cell(size(resblock));
    if strcmpi(outputQuantity_type,'vector')
        fprintf('\t>>   The quantity %s is a vector.\n',outputQuantity_name);
        for k = 1:numel(resblock)
            out{k} = textscan(resblock{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
            out{k} = horzcat(out{k}{:});
        end
    elseif strcmpi(outputQuantity_type,'scalar')
        fprintf('\t>>   The quantity %s is a scalar.\n',outputQuantity_name);
        for k = 1:numel(resblock)
            out{k} = textscan(resblock{k},'%f %f','delimiter',' ','MultipleDelimsAsOne', 1);
            out{k} = horzcat(out{k}{:});
        end
    else
        fprintf('\t>>   Error: The type of the quantity %s is not supported yet.\n',outputQuantity_name);
    end
    %
    if strcmp(outputQuantity_location,'OnNodes')
        outputQuantity.dataOnNodes(:,:,nbquantityperstep) = cell2mat(out);
    elseif strcmp(outputQuantity_location,'OnGaussPoints')
        outputQuantity.dataOnGaussPoints(:,:,nbquantityperstep) = cell2mat(out);
    else
        fprintf('\t>> Error: Location not supported.\n');
        return;
    end
    outputQuantity_step_old=outputQuantity_step;
end
%
%% Write tecplot file for the last step or if there is only 1 step
fprintf('\t>>   Flush the last step in a Tecplot file....\n');
writeData2plt(basename,nameFileMshInfo,outputQuantity_step_old, ...
    outputQuantity,headersMsh_Nnode,num_nodes, ...
    num_elements,msh,nbquantityperstep,ntotstep)
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData2plt(basename,nameFileMshInfo,outputQuantity_step_old, ...
    outputQuantity,headersMsh_Nnode,num_nodes, ...
    num_elements,msh,nbquantityperstep,ntotstep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write data into tecplot files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization and setting of variables
fieldwidth=20;
sigfigures=10;
numformat = strcat('%',num2str(fieldwidth),'.',num2str(sigfigures),'e');
% use a fixed length format for the step in the tecplot file
stepformat=strcat('%0',num2str(numel(num2str(ntotstep))),'d');
breaklineCharNum = 4000;
dotlocation = find(nameFileMshInfo =='.');
%
%% create the name of the tecplot file for the step and open it
nameFileMeshInfo_tecplot = strcat(nameFileMshInfo(1:dotlocation(end)-1),'_step_',num2str(outputQuantity_step_old,stepformat),'.dat');
fileID = fopen(nameFileMeshInfo_tecplot,'w');
%
%% create the headers of the tecplot file
list_tecplot_var = strjoin(outputQuantity.varname,' ');
list_tecplot_varlocation = strjoin(outputQuantity.varlocation,' ');
if headersMsh_Nnode == 3
    header_tecplot=strcat('TITLE = "',basename,'"\n VARIABLES = "X" "Y" "Z"'," ", ...
        list_tecplot_var,"\n ", ...
        'ZONE T="STEP %d"\n N=%d, E=%d, DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE, VARLOCATION = (NODAL NODAL NODAL '," ", ...
        list_tecplot_varlocation,')\n');
elseif headersMsh_Nnode == 4
    header_tecplot=strcat('TITLE = "',basename,'"\n VARIABLES = "X" "Y" "Z"'," ", ...
        list_tecplot_var,"\n ", ...
        'ZONE T="STEP %d"\n N=%d, E=%d, DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL, VARLOCATION = (NODAL NODAL NODAL '," ", ...
        list_tecplot_varlocation,')\n');
end
fprintf(fileID,header_tecplot,outputQuantity_step_old,num_nodes,num_elements);
%
%% coordinate (x,y,z) loop
for i=1:3
    % node loop
    for j=1:num_nodes
        fprintf(fileID, numformat, msh.nodes(j,i));
        if mod((j*fieldwidth), breaklineCharNum) == 0
            fprintf(fileID, '\n'); % Start a new line every 4000 characters
        end
    end
    fprintf(fileID,'\n');
end
%
%% output quantity loop
for iquantity=1:nbquantityperstep
    if strcmpi(outputQuantity.type(iquantity),'vector')
        % vector quantity
        if strcmp(outputQuantity.location(iquantity),'OnNodes')
            orderquantityOnNodes = zeros(length(outputQuantity.dataOnNodes(:,1,iquantity)),4);
            %orderquantityOnNodes = sortrows(outputQuantity);  % NOTE the multiplication by a scaling factor!!!
            orderquantityOnNodes = [outputQuantity.dataOnNodes(:,2:4,iquantity), outputQuantity.dataOnNodes(:,1,iquantity)]; % swap node ID to the last column
            for i=1:3
                % node loop
                for j=1:num_nodes
                    fprintf(fileID, numformat, orderquantityOnNodes(j,i));
                    if mod((j*fieldwidth), breaklineCharNum) == 0
                        fprintf(fileID, '\n'); % Start a new line every 4000 characters
                    end
                end
                fprintf(fileID,'\n');
            end
        elseif strcmp(outputQuantity.location(iquantity),'OnGaussPoints')
            orderquantityOnGaussPoints = zeros(length(outputQuantity.dataOnGaussPoints(:,1,iquantity)),4);
            %orderquantityOnGaussPoints = sortrows(outputQuantity);  % NOTE the multiplication by a scaling factor!!!
            orderquantityOnGaussPoints = [outputQuantity.dataOnGaussPoints(:,2:4,iquantity), outputQuantity.dataOnGaussPoints(:,1,iquantity)]; % swap node ID to the last column
            for i=1:3
                % element loop
                for j=1:num_elements
                    fprintf(fileID, numformat, orderquantityOnGaussPoints(j,i));
                    if mod((j*fieldwidth), breaklineCharNum) == 0
                        fprintf(fileID, '\n'); % Start a new line every 4000 characters
                    end
                end
                fprintf(fileID,'\n');
            end
        else
            fprintf('\t>> Error: Location not supported.\n');
            return
        end % end test on variable location
    else
        % scalar quantity
        if strcmp(outputQuantity.location(iquantity),'OnNodes')
            orderquantityOnNodes = zeros(length(outputQuantity.dataOnNodes(:,1,iquantity)),2);
            %orderquantityOnNodes = sortrows(outputQuantity);  % NOTE the multiplication by a scaling factor!!!
            orderquantityOnNodes = [outputQuantity.dataOnNodes(:,2,iquantity), outputQuantity.dataOnNodes(:,1,iquantity)]; % swap node ID to the last column
            % node loop
            for j=1:num_nodes
                fprintf(fileID, numformat, orderquantityOnNodes(j,1));
                if mod((j*fieldwidth), breaklineCharNum) == 0
                    fprintf(fileID, '\n'); % Start a new line every 4000 characters
                end
            end
            fprintf(fileID,'\n');
        elseif strcmp(outputQuantity.location(iquantity),'OnGaussPoints')
            orderquantityOnGaussPoints = zeros(length(outputQuantity.dataOnGaussPoints(:,1,iquantity)),2);
            %orderquantityOnGaussPoints = sortrows(outputQuantity);  % NOTE the multiplication by a scaling factor!!!
            orderquantityOnGaussPoints = [outputQuantity.dataOnGaussPoints(:,2,iquantity), outputQuantity.dataOnGaussPoints(:,1,iquantity)]; % swap node ID to the last column
            % element loop
            for j=1:num_elements
                fprintf(fileID, numformat, orderquantityOnGaussPoints(j,1));
                if mod((j*fieldwidth), breaklineCharNum) == 0
                    fprintf(fileID, '\n'); % Start a new line every 4000 characters
                end
            end
            fprintf(fileID,'\n');
        else
            fprintf('\t>> Error: Location not supported.\n');
            return
        end % end test on variable location
    end % end of test on outputQuantity.type
end % end of output quantity loop
%
%% connectivity loop
if headersMsh_Nnode == 3
    for i=1:num_elements
        for j=1:3
            %fprintf(fileID, '%10d',msh.elements(i,j));
            fprintf(fileID, '%10d',msh.nodeID2nodeLine(msh.elements(i,j)));
        end
        fprintf(fileID,'\n');
    end
elseif headersMsh_Nnode == 4
    for i=1:num_elements
        for j=1:4
            %fprintf(fileID, '%10d',msh.elements(i,j));
            fprintf(fileID, '%10d',msh.nodeID2nodeLine(msh.elements(i,j)));
        end
        fprintf(fileID,'\n');
    end
else
    fprintf('\t>> Error: ElementType not supported.\n');
    return
end % end of connectivity loop
%
%% close the file
fclose(fileID);
%
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
