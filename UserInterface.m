function varargout = UserInterface(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UserInterface_OpeningFcn, ...
                   'gui_OutputFcn',  @UserInterface_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


function UserInterface_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

movegui(gcf,'center')

if(~isfield(handles, 'network'))
    set(handles.btnBrowse, 'Visible', 'off');
    set(handles.btnProcess, 'Visible', 'off');
    set(handles.resetBtn, 'Visible', 'off');
    set(handles.backBtn, 'Visible', 'off');
    set(handles.txtUrl, 'Visible', 'off');     
    set(handles.text36, 'Visible', 'off');    
    set(handles.text25, 'Visible', 'off');
    set(handles.text50, 'Visible', 'off');
    set(handles.text75, 'Visible', 'off');
    set(handles.text100, 'Visible', 'off');
    
    set(handles.origIm, 'Visible', 'off');
    set(handles.uipanel2, 'Visible', 'off');      
end
guidata(hObject, handles);


function varargout = UserInterface_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function txtUrl_Callback(hObject, eventdata, handles)


function txtUrl_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function btnBrowse_Callback(hObject, eventdata, handles)
handles.output = hObject;
if(isfield(handles, 'network'))
    [fn pn]=uigetfile('*jpg', 'select .jpg file');
    complete = strcat(pn, fn);
    if(complete)
        set(handles.txtUrl,'string', complete);
        
        handles.no = complete(length(complete)-6:length(complete)-4);
        %msgbox(num2str(complete(28:30)));
        handles.I = imread(complete);
        axes(handles.origIm);
        imshow(handles.I);
    end
end
guidata(hObject,handles);


function btnProcess_Callback(hObject, eventdata, handles)
handles.output = hObject;
if(isfield(handles, 'I'))
    img = handles.I;
    network = handles.network;
    
    if(isempty(get(handles.origIm, 'children')))
        imshow(img, 'Parent', handles.origIm);
    end

    datarsz = imresize(img, [512 512]);
    sizes = size(datarsz);

    datadbl = im2double (datarsz);
    datagry = rgb2gray (datadbl);
    handles.gry = datagry;

    datahist = histeq(datagry);
    datamed = medfilt2(datahist);
    d1ca = imadjust(datamed); 
    d1nF = wiener2(d1ca);
    imshow(d1nF, 'Parent', handles.preImg);

    b = 255;
    trans_thresh = mrdivide(0:1:b,b);                                   

    for j=1:(b+1)
        value_trans1 = trans_thresh(1,j);
        result4_trans1(j,1) = bweuler(im2bw(datahist,value_trans1),4);
    end

    max_result4_trans1 = max(result4_trans1);
    min_result4_trans1 = min(result4_trans1);
    m = 0;
    n = 0;
    max_sum_result4_trans1 = 0;
    min_sum_result4_trans1 = 0;

    for j = 1:(b+1)
        if result4_trans1(j,1) == max_result4_trans1
            m = m+1;
            max_sum_result4_trans1 = max_sum_result4_trans1 + j;
        elseif result4_trans1(j,1) == min_result4_trans1
            n = n+1;
            min_sum_result4_trans1 = min_sum_result4_trans1 + j;    
        else
        end
    end

    threshold_I1 = ((max_sum_result4_trans1)+(min_sum_result4_trans1))/(m+n);
    dataseg =  im2bw(datahist,threshold_I1/(b+1));
    
    datangv = imcomplement(dataseg);
    %test = edges(datamed, 'sobel');

    datafil = imfill(datangv, 'holes');

    se = strel('disk',10);
    datarode = imerode(datafil, se);

    path = ['masks/mask' handles.no '.jpg'];
    %msgbox(path);
    maskS = imread(path);
    maskS = im2bw(maskS);

    [labeled, num] = bwlabel(maskS);

    temp = zeros(512,512);

    for n=1:num
        [i,j] = find(labeled==n,1,'first');
        labelS = bwselect(datarode,j,i,4);
        temp = temp | labelS;    
    end
   
    datadlt = imdilate(temp, se);

    se = strel('line',25,140); % 74 degrees determined by inspection
    bw2 = imclose(datadlt,se);
    se2 = strel('line',25,140+90);
    datasmt = imclose(bw2,se2);

    t = handles.gry;
    [m,n] = size(datasmt);
    idx = find(datasmt==1);
    hasil = zeros(m,n);
    hasil(idx) = t(idx);
    %hasil = uint8(hasil);

    imshow(dataseg, 'Parent', handles.seg1Img);
    imshow(temp, 'Parent', handles.binImg);
    imshow(hasil, 'Parent', handles.seg2Img);
    %%
    p = imhist(hasil);
    p = p./numel(hasil);
    L = length(p);
    [v, mu] = statmoments(p, 3);

    glcm0 = graycomatrix(hasil, 'Offset', [0 1]);
    stats0 = graycoprops(glcm0,{'contrast'});
    arrays0 = struct2array(stats0);

    glcm45 = graycomatrix(hasil, 'Offset', [-1 1]);
    stats45 = graycoprops(glcm45,{'contrast'});
    arrays45 = struct2array(stats45);

    glcm90 = graycomatrix(hasil, 'Offset', [-1 0]);
    stats90 = graycoprops(glcm90,{'contrast'});
    arrays90 = struct2array(stats90);

    glcm135 = graycomatrix(hasil, 'Offset', [-1 -1]);
    stats135 = graycoprops(glcm135,{'contrast'});
    arrays135 = struct2array(stats135);

    mcontrast = [arrays0 arrays45 arrays90 arrays135];
    contrast = mean(mcontrast);
    set(handles.txtContrast, 'String', contrast);

    %%
    eglcm0 = graycomatrix(hasil, 'Offset', [0 1]);
    estats0 = graycoprops(eglcm0,{'energy'});
    earrays0 = struct2array(estats0);

    eglcm45 = graycomatrix(hasil, 'Offset', [-1 1]);
    estats45 = graycoprops(eglcm45,{'energy'});
    earrays45 = struct2array(estats45);

    eglcm90 = graycomatrix(hasil, 'Offset', [-1 0]);
    estats90 = graycoprops(eglcm90,{'energy'});
    earrays90 = struct2array(estats90);

    eglcm135 = graycomatrix(hasil, 'Offset', [-1 -1]);
    estats135 = graycoprops(eglcm135,{'energy'});
    earrays135 = struct2array(estats135);

    menergy = [earrays0 earrays45 earrays90 earrays135];
    energy = mean(menergy);
    set(handles.txtEnergy, 'String', energy);

    %%
    avg = mean(hasil(:));
    t(1) = mu(1);
    set(handles.txtAVG, 'String', t(1)); 
    
    %%
    stdev = std(hasil(:));
    t(2) = mu(2).^0.5;
    set(handles.txtSD, 'String', t(2));

    %%
    vars = var(hasil(:));
    t(3) = 1-(1/(1+vars));        
    set(handles.txtSmoothness, 'String', t(3));        

    %%
    skew = skewness(hasil(:));
    t(4) = skew;
    set(handles.txtThirdMoment, 'String', t(4));        

    %%
    uniform = (1-(stdev/avg));
    t(5) = sum(p.^2);
    set(handles.txtUniformity, 'String', t(5));        

    
    %%
    entrop = entropy(hasil);
    t(6) = -sum(p.*(log2(p + eps)));
    set(handles.txtEntropy, 'String', t(6));



    %%
    new_inputs=[t(1) t(2) t(3) t(4) t(5) t(6) contrast energy]';

    namafile=get(handles.txtUrl,'String') 
    target = zeros(1,2);
    if(strfind(namafile, 'BGN'))
        target(1, 1) = 1;
        target(1, 2) = 0;
        targets = 0;
    elseif(strfind(namafile, 'MGN'))
        target(1, 1) = 0;
        target(1, 2) = 1;
        targets = 1;
    end
    new_targets=target'
    targets
    netww = handles.network;
    [result,c]=classification(netww, new_inputs, new_targets,targets);
    set(handles.txtResult, 'String',result);
   
else
    msgbox('Please choose file first');
end
guidata(hObject,handles);


function resetBtn_Callback(hObject, eventdata, handles)
cla(handles.preImg, 'reset');
cla(handles.seg2Img, 'reset');
cla(handles.binImg, 'reset');
cla(handles.seg1Img, 'reset');
cla(handles.origIm,'reset');
        
set(handles.txtAVG, 'String', '');
set(handles.txtSD, 'String', '');
set(handles.txtSmoothness, 'String', '');
set(handles.txtThirdMoment, 'String', '');
set(handles.txtUniformity, 'String', '');
set(handles.txtEntropy, 'String', '');
set(handles.txtResult, 'String','');  
set(handles.txtContrast, 'String','');  
set(handles.txtEnergy, 'String','');  


function backBtn_Callback(hObject, eventdata, handles)
close(handles.figure1)


% --- Executes on button press in btnTraining.
function btnTraining_Callback(hObject, eventdata, handles)
    imgPath = '.\training'; 
    %h = waitbar(40,'Training data...');
    
    files = dir(fullfile(imgPath, '*.jpg')); 
    nfiles = length(files);
    inputs = zeros(nfiles,8);
    target = zeros(nfiles,2);
    set(handles.text36, 'Visible', 'on');    
    fprintf('-------------------------\n');
    
    %tstart = tic
       
    for ii=1:nfiles
        if(ii == nfiles*0.25)
            set(handles.text36, 'Visible', 'off');            
            set(handles.text25, 'Visible', 'on');            
        elseif(ii == nfiles*0.50)
            set(handles.text25, 'Visible', 'off');
            set(handles.text50, 'Visible', 'on');
        elseif(ii == nfiles*0.75)
            set(handles.text75, 'Visible', 'on');
            set(handles.text50, 'Visible', 'off');            
        end
        drawnow();
    
        currentfilename = [imgPath '\' files(ii).name];
        currentimage = im2double(imread(currentfilename));
        images{ii} = currentimage;
        data = rgb2gray(images{ii});
        handles.gry = data;

        check = 0;
        if strfind(currentfilename, 'BGN') 
            check = check+1;
        elseif strfind(currentfilename, 'MGN') 
            check = check+1;
        end

        if check == 1
            fprintf('%d) Input file: %s\n', ii, currentfilename);
            %---------------
            %fitur ekstraksi
            %---------------
            datarsz = imresize(data, [512 512]);
            sizes = size(datarsz);

            datahist = histeq(datarsz);
            datamed = medfilt2(datahist);
            d1ca = imadjust(datamed); 
            d1nF = wiener2(d1ca);

            b = 255;
            trans_thresh = mrdivide(0:1:b,b);                                  
            for j=1:(b+1)
                value_trans1 = trans_thresh(1,j);
                result4_trans1(j,1) = bweuler(im2bw(datahist,value_trans1),4);
            end

            max_result4_trans1 = max(result4_trans1);
            min_result4_trans1 = min(result4_trans1);
            m = 0;
            n = 0;
            max_sum_result4_trans1 = 0;
            min_sum_result4_trans1 = 0;

            for j = 1:(b+1)
            if result4_trans1(j,1) == max_result4_trans1
                m = m+1;
                max_sum_result4_trans1 = max_sum_result4_trans1 + j;
            elseif result4_trans1(j,1) == min_result4_trans1
                n = n+1;
                min_sum_result4_trans1 = min_sum_result4_trans1 + j;    
            else
            end
            end

            threshold_I1 = ((max_sum_result4_trans1)+(min_sum_result4_trans1))/(m+n);
            dataseg =  im2bw(datahist,threshold_I1/(b+1));
            
            datangv = imcomplement(dataseg);
            %test = edges(datamed, 'sobel');

            datafil = imfill(datangv, 'holes');

            se = strel('disk',10);
            datarode = imerode(datafil, se);

            path = ['masks/mask' num2str(currentfilename(length(currentfilename)-6:length(currentfilename)-4)) '.jpg'];
           % msgbox(path);
            maskS = imread(path);
            maskS = im2bw(maskS);

            [labeled, num] = bwlabel(maskS);

            temp = zeros(512,512);

            for n=1:num
                [i,j] = find(labeled==n,1,'first');
                labelS = bwselect(datarode,j,i,4);
                temp = temp | labelS;    
            end

            datadlt = imdilate(temp, se);

            se = strel('line',25,140); % 74 degrees determined by inspection
            bw2 = imclose(datadlt,se);
            se2 = strel('line',25,140+90);
            datasmt = imclose(bw2,se2);

            t = handles.gry;
            [m,n] = size(datasmt);
            idx = find(datasmt==1);
            hasil = zeros(m,n);
            hasil(idx) = t(idx);
            %hasil = uint8(hasil);

            
            p = imhist(hasil);
            p = p./numel(hasil);
            L = length(p);
            [v, mu] = statmoments(p, 3);

            
            %%
            t(1) = mu(1);

            %%
            %stdev = std(hasil(:));
            t(2) = mu(2).^0.5;
            
            %%
            vars = var(hasil(:));
            t(3) = 1-(1/(1+vars));        
            
            
            %%
            t(4) = skewness(hasil(:));
            
            %%
            %uniform = (1-(stdev/avg));
            t(5) = sum(p.^2);
            
            %%
            %entrop = entropy(hasil);
            t(6) = -sum(p.*(log2(p + eps)));            
            
            %%
            glcm0 = graycomatrix(hasil, 'Offset', [0 1]);
            stats0 = graycoprops(glcm0,{'contrast'});
            arrays0 = struct2array(stats0);

            glcm45 = graycomatrix(hasil, 'Offset', [-1 1]);
            stats45 = graycoprops(glcm45,{'contrast'});
            arrays45 = struct2array(stats45);

            glcm90 = graycomatrix(hasil, 'Offset', [-1 0]);
            stats90 = graycoprops(glcm90,{'contrast'});
            arrays90 = struct2array(stats90);

            glcm135 = graycomatrix(hasil, 'Offset', [-1 -1]);
            stats135 = graycoprops(glcm135,{'contrast'});
            arrays135 = struct2array(stats135);

            mcontrast = [arrays0 arrays45 arrays90 arrays135];
            vcontrast = mean(mcontrast);


            %%
            eglcm0 = graycomatrix(hasil, 'Offset', [0 1]);
            estats0 = graycoprops(eglcm0,{'energy'});
            earrays0 = struct2array(estats0);

            eglcm45 = graycomatrix(hasil, 'Offset', [-1 1]);
            estats45 = graycoprops(eglcm45,{'energy'});
            earrays45 = struct2array(estats45);

            eglcm90 = graycomatrix(hasil, 'Offset', [-1 0]);
            estats90 = graycoprops(eglcm90,{'energy'});
            earrays90 = struct2array(estats90);

            eglcm135 = graycomatrix(hasil, 'Offset', [-1 -1]);
            estats135 = graycoprops(eglcm135,{'energy'});
            earrays135 = struct2array(estats135);

            menergy = [earrays0 earrays45 earrays90 earrays135];
            venergy = mean(menergy);
            
            %%
            inputs(ii, 1) = t(1);
            inputs(ii, 2) = t(2);
            inputs(ii, 3) = t(3);
            inputs(ii, 4) = t(4); 
            inputs(ii, 5) = t(5);
            inputs(ii, 6) = t(6);
            inputs(ii, 7) = vcontrast;  
            inputs(ii, 8) = venergy;  

            if(strfind(currentfilename, 'BGN'))
                target(ii, 1) = 1;
                target(ii, 2) = 0;
            elseif(strfind(currentfilename, 'MGN'))
                target(ii, 1) = 0;
                target(ii, 2) = 1;
            end
        end
    end
    
    
    %tend = toc(tstart)
    
    times = ['Training data completed'];
    set(handles.text38, 'String', times);   
    
    set(handles.text75, 'Visible', 'off');    
    set(handles.text100, 'Visible', 'on');
    drawnow();
    %msgbox(times, 'Completed msgbox')

    %-------------
    %training
    %-------------
    inputs = inputs'
    targets = target'
    hiddenLayerSize = 6;
    net = patternnet(hiddenLayerSize);
    net.divideParam.trainRatio = 100/100;
    net.trainParam.showWindow = false;
    net.trainParam.showCommandLine = false; 
    [net,tr] = train(net,inputs,targets);   
    handles.network = net;    
    
    
    set(handles.btnBrowse, 'Visible', 'on');
    set(handles.btnProcess, 'Visible', 'on');
    set(handles.resetBtn, 'Visible', 'on');

    set(handles.backBtn, 'Visible', 'on');
    set(handles.txtUrl, 'Visible', 'on');
    
    set(handles.origIm, 'Visible', 'on');
    set(handles.uipanel2, 'Visible', 'on');     
    set(handles.btnTraining, 'String', 'Trained', 'enable', 'off');
    
    
guidata(hObject,handles);


% --- Executes on button press in btnTraining.
