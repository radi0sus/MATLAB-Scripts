% Analyzes orbitals from ORCA output files.
% ORCA: https://orcaforum.kofo.mpg.de
% Contributions of certain atoms or elements to orbitals are investigated.
% At the moment only 'LOEWDIN REDUCED ORBITAL POPULATIONS PER MO' are supported.
% 
% Use the '!LargePrint' option of ORCA.
%
% 1. Check the options.
% 2. Run the script with F5. 
% 3. Open an ORCA output file.
% 4. Output will be written in 'orb-analysis.txt' in the same folder
%    - You need write permissions!
%    - Previous 'orb-analysis.txt' file will be overwritten!
%
% 
% options:
%
% orbitals_to_analyze:
% --------------------
%
%  The orbitals_to_analyze instruction recognizes only the keywords 
%  below or integers!
%
% 'HOMO+-number' = analyze +- number of orbitals from HOMO
%  example 1: orbitals_to_analyze='HOMO+-10'
%             If the orbital number of the HOMO is 50, the script analyzes
%             every orbital in the range from 40 to 60.
%  example 2: orbitals_to_analyze='HOMO+-0'
%             Only the HOMO will be analyzed.
%
% Number = analyze orbital with the given Number
% example: orbitals_to_analyze=50
%          The orbital number 50 will be analyzed
%
% Number:Number = analyze orbitals within the given range
% example: orbitals_to_analyze=0:50
%          Every orbital in the range from 0 to 50 will be analyzed.
%
% 'all' = analyze all orbitals
%  example: orbitals_to_analyze='all'
%           All orbitals will be analyzed.
%
% threshold:
% --------------------
% 
% 'threshold' = print orbitals above the given threshold
%  example 1: threshold = 5.2
%             Every atom contribution to a specific orbital will be
%             printed in the output file, if the threshold of the 
%             contribution is greater than 5.2.
%            
%  example 2: threshold = -1
%             Every atom contribution to a specific orbital listed in 
%             the ORCA output file will be printed.
%             
% 
% beta_orbitals:
% --------------------   
%
% 'beta_orbitals' = include beta orbitals (spin down) in the output
%  beta_orbitals = 1 % include beta orbitals in the output file
%  beta_orbitals = 0 % do not include beta orbitals in the output file
%
%

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% options
orbitals_to_analyze = 'HOMO+-10';
%orbitals_to_analyze = 1:10;
%orbitals_to_analyze = 'all';
threshold = 5;
beta_orbitals = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end options

% variables
orb_en_occ = {};
A = [];
O = [];
B = [];
k = 1;
l = 0;
spin_down_detected = 0;
Thcomp=[];
Thecomp=[];
Thbcomp=[];
Thbecomp=[];
hashmark=[];
% end variables

% open the output file
[orca_out_file, path] = uigetfile('.out','Select an ORCA output file');
if isequal(orca_out_file,0)
    error('Select an ORCA output file!');
    return;
else
    orca_out_file_ID=fopen(fullfile(path,orca_out_file),'r');
end

% look for string 'ORBITAL ENERGIES' in orca.out
% determine the total number of orbitals and the number of the HOMO
while ~feof(orca_out_file_ID)
    line=fgetl(orca_out_file_ID);
    is_energies=strfind(line,'ORBITAL ENERGIES');
    if is_energies
        % skip next 4 lines
        for nlines_2_skip = 1:4
            line=fgetl(orca_out_file_ID);
        end
        while ~isempty(line)
            % if there is no empty line but a '---' line, e.g. after restr. orbitals
            if strncmp(line,'-',1)
                break
            end
            %
            orb_en_occ(end+1,:) = textscan(line,'%f %f %*[^\n]');
            % determine HOMO number
            % only the spin up HOMO is considered
            if orb_en_occ{end,2} > 0
                homo_nr=orb_en_occ{end,1};
            end
            line = fgetl(orca_out_file_ID);
        end
        % determine total number of orbitals
        tot_nr_orb=length(orb_en_occ);
        break
    end
end

% look for string 'LOEWDIN REDUCED ORBITAL POPULATIONS PER MO' in orca.out
% in case of two (or more), take data from the last occurrence
while ~feof(orca_out_file_ID)
    line=fgetl(orca_out_file_ID);
    is_loewdin=strcmp(line,'LOEWDIN REDUCED ORBITAL POPULATIONS PER MO');
    if is_loewdin
        is_loewdin;
        % save the position of the last occurrence
        loewdin_pos=ftell(orca_out_file_ID);
    end
end

% move to the position
fseek(orca_out_file_ID,loewdin_pos,'bof');

% read all orbitals into struct
if isempty(is_loewdin)
    error('LOEWDIN REDUCED ORBITAL POPULATIONS PER MO not found!') % error message > end prgm
else
    while ~feof(orca_out_file_ID)           % loop to end of file
        while ~feof(orca_out_file_ID)       % loop to end of file
            line = fgetl(orca_out_file_ID); % read line by line
            % leave loop in case of an empty line
            if isempty(line)
                break
            end
            % check if spin down orbitals present
            if strfind(line,'SPIN DOWN')
                spin_down_detected = 1;
            end
            % remove some lines
            if strncmp(line,'THRESHOLD',1) || strncmp(line,'---------',1) || strncmp(line,'SPIN UP',1)
                continue
            end
            % remove ------ ------ lines
            if strfind(line,'--------')
                line = fgetl(orca_out_file_ID);
            end
            % split the line
            split_line = strsplit(line);
            % delete empty cells
            split_line = split_line(~cellfun(@isempty, split_line));
            % first 3 lines contain orb-nr, orb-en & orb-occ
            if k <= 3
                O = [O; split_line];
                k = k + 1;
            else
                % next lines contain populations per MO
                A = [A; split_line];
            end
        end
        % leave loop at the very last line -> else error
        if isempty(O)
            break
        end
        % struct: orbital(orb-nr).orb-en, orbital(orb-nr).occ,
        % orbital(orb-nr).T = table with population details
        % orbitals start with 1 (ORCA orb nr. 0 = 1)
        % spin up (alpha) or closed shell orbitals
        if ~spin_down_detected
            for i = str2double(O(1,1))+1 : str2double(O(1,1)) + length(O(1,:))
                l = l + 1;
                % extract the first letter of the orbital
                for n = 1:length(A(:,3))
                    B = [B;string(A{n,3}(1))];
                end
                orbital_a(i).en = str2double(O(2,l));
                orbital_a(i).occ = str2double(O(3,l));
                orbital_a(i).T = table(str2double(A(:,1)),string(A(:,2)),string(A(:,3)),B,str2double(A(:,l+3)),'VariableNames',{'aNum','aName','aOrb','aOrbs','aX'});
                B = [];
            end
            % reset all arrays and variables
            A = [];
            O = [];
            k = 1;
            l = 0;
        end
        % struct: orbital(orb-nr).orb-en, orbital(orb-nr).occ,
        % orbital(orb-nr).T = table with population details
        % orbitals start with 1 (ORCA orb nr. 0 = 1)
        % spin down (beta) orbitals
        if spin_down_detected
            for i = str2double(O(1,1))+1 : str2double(O(1,1)) + length(O(1,:))
                l = l+1;
                % extract the first letter of the orbital
                for n = 1:length(A(:,3))
                    B = [B;string(A{n,3}(1))];
                end
                orbital_b(i).en = str2double(O(2,l));
                orbital_b(i).occ = str2double(O(3,l));
                orbital_b(i).T = table(str2double(A(:,1)),string(A(:,2)),A(:,3),B,str2double(A(:,l+3)),'VariableNames',{'aNum','aName','aOrb','aOrbs','aX'});
                B = [];
            end
            % reset all arrays and variables
            A = [];
            O = [];
            k = 1;
            l = 0;
        end
    end
end
fclose(orca_out_file_ID);
% end reading orbitals


% analyze orbitals

% check options, which orbitals should be printed
if strncmp(orbitals_to_analyze,'HOMO+-',6)
    orbitals_to_analyze = str2num(orbitals_to_analyze(7:end));
    loop_start = homo_nr - orbitals_to_analyze;
    loop_end = homo_nr + orbitals_to_analyze;
elseif isnumeric(orbitals_to_analyze)
    loop_start = orbitals_to_analyze(1);
    loop_end = orbitals_to_analyze(end);
elseif strncmp(orbitals_to_analyze,'all',3)
    loop_start = 0;
    loop_end = tot_nr_orb-1;
else
    error('Keyword not recognized!');
    return
end

% save analysis
% open file
analysis_out_file=[path 'orb-analysis.txt'];
fID = fopen(analysis_out_file,'w');

fprintf(fID,'-------------------------------------------------------------------------------\n');
fprintf(fID,'orbital analysis of %s\n',orca_out_file);

if spin_down_detected
    fprintf(fID,'total no. of orbitals: %d\n',(tot_nr_orb-1)*2);
else
    fprintf(fID,'total no. of orbitals: %d\n',tot_nr_orb-1);
end
fprintf(fID,'HOMO no.: %d\n',homo_nr);
fprintf(fID,'orbitals to analyze: %d...%d\n',loop_start,loop_end);
if spin_down_detected
    fprintf(fID,'beta orbitals detected: yes\n');
else
    fprintf(fID,'beta orbitals detected: no\n');
end
if beta_orbitals && spin_down_detected
    fprintf(fID,'beta orbitals included: yes\n');
else
    fprintf(fID,'beta orbitals included: no\n');
end
fprintf(fID,'threshold for printing orbitals (%%): %d\n',threshold);
fprintf(fID,'-------------------------------------------------------------------------------\n');
fprintf(fID,' \n');

for i = loop_start:loop_end
    % check if range of orbitals is ok
    if loop_start < 0 || i < 0 || i >  tot_nr_orb-1 || loop_end > tot_nr_orb-1
        error('Values exceed range of orbitals!');
        return
    end
    % end check
    fprintf(fID,'*******************************************************************************\n');
    fprintf(fID,'orbital no: %d, energy: %f, occ: %d\n',i,orbital_a(i+1).en,orbital_a(i+1).occ);
    fprintf(fID,'*******************************************************************************\n');
    Tecomp=varfun(@sum,orbital_a(i+1).T,'InputVariables','aX','GroupingVariables',{'aName'});
    Tecomp=Tecomp(Tecomp{:,3} > threshold,:);
    Tecomp.Properties.VariableNames={'Element', 'GroupCount','Contribution'};
    Tacomp=varfun(@sum,orbital_a(i+1).T,'InputVariables','aX','GroupingVariables',{'aNum','aName'});
    Tacomp=Tacomp(Tacomp{:,4} > threshold,:);
    Tacomp.Properties.VariableNames={'AtomNo', 'Element', 'GroupCount','Contribution'};
    Tocomp=varfun(@sum,orbital_a(i+1).T,'InputVariables','aX','GroupingVariables',{'aNum','aName','aOrbs'});
    Tocomp=Tocomp(Tocomp{:,5} > threshold,:);
    Tocomp.Properties.VariableNames={'AtomNo', 'Element', 'Orbital', 'GroupCount','Contribution'};
    Tgcomp=orbital_a(i+1).T;
    Tgcomp=removevars(Tgcomp,'aOrbs');
    Tgcomp=Tgcomp(Tgcomp{:,4} > threshold,:);
    Tgcomp.Properties.VariableNames={'AtomNo', 'Element', 'Orbital','Contribution'};
    %%% for the summary
    Thmcomp=varfun(@sum,orbital_a(i+1).T,'InputVariables','aX','GroupingVariables',{'aName'});
    [numr, numc] = size(Thmcomp);
    orbNum = [repmat(i,1,numr)].';
    orbOcc = [repmat(orbital_a(i+1).occ,1,numr)].';
    orbEn = [repmat(orbital_a(i+1).en,1,numr)].';
    Thmcomp = addvars(Thmcomp,orbNum,orbOcc,orbEn,'Before','aName');
    for k = 1:numr
        hashmark = [string(hashmark);string(repmat('*',1,round(Thmcomp.sum_aX(k)/10)))];
    end
    Thmcomp = addvars(Thmcomp,hashmark,'After','sum_aX');
    Thmcomp = removevars(Thmcomp,'GroupCount');
    Thmcomp = Thmcomp(Thmcomp{:,5} > threshold,:);
    Thcomp = [Thcomp;{"","","","","",""};Thmcomp];
    Thecomp = [Thecomp;Thmcomp];
    hashmark=[];
    %%%
    % this tables will be finally printed in the output file
    a=[Tecomp.Properties.VariableNames;{'-------','----------','------------'};Tecomp{:,:}];
    fprintf(fID,' \n');
    fprintf(fID,'%+6s %+10s %+12s\n',a.');
    fprintf(fID,' \n');
    a=[Tacomp.Properties.VariableNames;{'------','-------','----------','------------'};Tacomp{:,:}];
    fprintf(fID,' \n');
    fprintf(fID,'%+5s %+6s %+10s %+12s\n',a.');
    fprintf(fID,' \n');
    a=[Tocomp.Properties.VariableNames;{'------','-------','-------','----------','------------'};Tocomp{:,:}];
    fprintf(fID,' \n');
    fprintf(fID,'%+5s %+6s %+7s %+10s %+12s\n',a.');
    fprintf(fID,' \n');
    a=[Tgcomp.Properties.VariableNames;{'------','-------','-------','------------'};Tgcomp{:,:}];
    fprintf(fID,' \n');
    fprintf(fID,'%+5s %+6s %+7s %+12s\n',a.');
    fprintf(fID,' \n');
end

if beta_orbitals && spin_down_detected
    for i = loop_start:loop_end
        % check if range of orbitals is ok
        if i < 0 || i >  tot_nr_orb
            error('Values exceed range of orbitals!');
            return
        end
        % end check
        fprintf(fID,'*******************************************************************************\n');
        fprintf(fID,'orbital no: %d (beta), energy: %f, occ: %d\n',i,orbital_b(i+1).en,orbital_b(i+1).occ);
        fprintf(fID,'*******************************************************************************\n');
        Tecomp=varfun(@sum,orbital_b(i+1).T,'InputVariables','aX','GroupingVariables',{'aName'});
        Tecomp=Tecomp(Tecomp{:,3} > threshold,:);
        Tecomp.Properties.VariableNames={'Element', 'GroupCount','Contribution'};
        Tacomp=varfun(@sum,orbital_b(i+1).T,'InputVariables','aX','GroupingVariables',{'aNum','aName'});
        Tacomp=Tacomp(Tacomp{:,4} > threshold,:);
        Tacomp.Properties.VariableNames={'AtomNo', 'Element', 'GroupCount','Contribution'};
        Tocomp=varfun(@sum,orbital_b(i+1).T,'InputVariables','aX','GroupingVariables',{'aNum','aName','aOrbs'});
        Tocomp=Tocomp(Tocomp{:,5} > threshold,:);
        Tocomp.Properties.VariableNames={'AtomNo', 'Element', 'Orbital', 'GroupCount','Contribution'};
        Tgcomp=orbital_b(i+1).T;
        Tgcomp=removevars(Tgcomp,'aOrbs');
        Tgcomp=Tgcomp(Tgcomp{:,4} > threshold,:);
        Tgcomp.Properties.VariableNames={'AtomNo', 'Element', 'Orbital','Contribution'};
        %%% for the summary
        Thmcomp=varfun(@sum,orbital_b(i+1).T,'InputVariables','aX','GroupingVariables',{'aName'});
        [numr, numc] = size(Thmcomp);
        orbNum = [repmat(i,1,numr)].';
        orbOcc = [repmat(orbital_b(i+1).occ,1,numr)].';
        orbEn = [repmat(orbital_b(i+1).en,1,numr)].';
        Thmcomp = addvars(Thmcomp,orbNum,orbOcc,orbEn,'Before','aName');
        for k = 1:numr
            hashmark = [string(hashmark);string(repmat('*',1,round(Thmcomp.sum_aX(k)/10)))];
        end
        Thmcomp = addvars(Thmcomp,hashmark,'After','sum_aX');
        Thmcomp = removevars(Thmcomp,'GroupCount');
        Thmcomp = Thmcomp(Thmcomp{:,5} > threshold,:);
        Thbcomp = [Thbcomp;{"","","","","",""};Thmcomp];
        Thbecomp = [Thbecomp;Thmcomp];
        hashmark=[];
        %%%
        % this tables will be finally printed in the output file
        a=[Tecomp.Properties.VariableNames;{'-------','----------','------------'};Tecomp{:,:}];
        fprintf(fID,' \n');
        fprintf(fID,'%+6s %+10s %+12s\n',a.');
        fprintf(fID,' \n');
        a=[Tacomp.Properties.VariableNames;{'------','-------','----------','------------'};Tacomp{:,:}];
        fprintf(fID,' \n');
        fprintf(fID,'%+5s %+6s %+10s %+12s\n',a.');
        fprintf(fID,' \n');
        a=[Tocomp.Properties.VariableNames;{'------','-------','-------','----------','------------'};Tocomp{:,:}];
        fprintf(fID,' \n');
        fprintf(fID,'%+5s %+6s %+7s %+10s %+12s\n',a.');
        fprintf(fID,' \n');
        a=[Tgcomp.Properties.VariableNames;{'------','-------','-------','------------'};Tgcomp{:,:}];
        fprintf(fID,' \n');
        fprintf(fID,'%+5s %+6s %+7s %+12s\n',a.');
        fprintf(fID,' \n');
    end
end
%%%
Thcomp.Properties.VariableNames={'No', 'Occ', 'Orbital_Energy', 'Element', 'Cntrb','Cntrb_X'};
Thcomp = movevars(Thcomp,{'Orbital_Energy'},'After','Cntrb_X');
fprintf(fID,'*******************************************************************************\n');
fprintf(fID,'summary of orbital contributions\n');
fprintf(fID,'*******************************************************************************\n');
a=[Thcomp.Properties.VariableNames;{'--','---','-------','------','----------','--------------'};Thcomp{:,:}];
fprintf(fID,' \n');
fprintf(fID,'%+5s %+3s %+7s %+6s %-10s %+14s\n',a.');
fprintf(fID,' \n');
fprintf(fID,'*******************************************************************************\n');
fprintf(fID,'summary of element contributions to orbitals\n');
fprintf(fID,'*******************************************************************************\n');
Thecomp.Properties.VariableNames={'No', 'Occ', 'Orbital_Energy', 'Element', 'Cntrb','Cntrb_X'};
Thecomp = sortrows(Thecomp,[2 4 5],{'descend' 'ascend' 'ascend'});
Thecomp = movevars(Thecomp,{'No','Occ','Orbital_Energy'},'After','Cntrb_X');
a=[Thecomp.Properties.VariableNames;{'-------','------','----------','--','---','--------------'};Thecomp{:,:}];
fprintf(fID,' \n');
fprintf(fID,'%+7s %+6s %-10s %+4s %+3s %+14s\n',a.');
fprintf(fID,' \n');

if beta_orbitals && spin_down_detected
    Thbcomp.Properties.VariableNames={'No', 'Occ', 'Orbital_Energy', 'Element', 'Cntrb','Cntrb_X'};
    Thbcomp = movevars(Thcomp,{'Orbital_Energy'},'After','Cntrb_X');
    fprintf(fID,'*******************************************************************************\n');
    fprintf(fID,'summary of orbital contributions (beta)\n');
    fprintf(fID,'*******************************************************************************\n');
    a=[Thbcomp.Properties.VariableNames;{'--','---','-------','------','----------','--------------'};Thbcomp{:,:}];
    fprintf(fID,' \n');
    fprintf(fID,'%+5s %+3s %+7s %+6s %-10s %+14s\n',a.');
    fprintf(fID,' \n');
    fprintf(fID,'*******************************************************************************\n');
    fprintf(fID,'summary of element contributions to orbitals (beta)\n');
    fprintf(fID,'*******************************************************************************\n');
    Thbecomp.Properties.VariableNames={'No', 'Occ', 'Orbital_Energy', 'Element', 'Cntrb','Cntrb_X'};
    Thbecomp = sortrows(Thbecomp,[2 4 5],{'descend' 'ascend' 'ascend'});
    Thbecomp = movevars(Thbecomp,{'No','Occ','Orbital_Energy'},'After','Cntrb_X');
    a=[Thbecomp.Properties.VariableNames;{'-------','------','----------','--','---','--------------'};Thbecomp{:,:}];
    fprintf(fID,' \n');
    fprintf(fID,'%+7s %+6s %-10s %+4s %+3s %+14s\n',a.');
    fprintf(fID,' \n');
end

%close file
fclose(fID);
