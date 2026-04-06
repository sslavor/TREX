
function TREX_UI

    % Main window
    fig = uifigure('Name','TREX – Turbine Rotor EXtractor',...
                   'Position',[100 100 1000 600],...
                   'Color',[0 0 0]);   % Dark background
    
    %% ===== RIGHT SIDE (Preview Image) =====
    
    previewPanel = uipanel(fig,...
        'Title','',...
        'FontSize',14,...
        'Position',[350 50 600 800],...
        'BackgroundColor',[0 0 0],...
        'HighlightColor',[0 0 0],...
        'BorderColor',[0 0 0]);

    
    img = uiimage(previewPanel,...
        'Position',[0 0 previewPanel.Position(3) previewPanel.Position(4)]);
    img.ImageSource = 'TREX.png';
    img.ScaleMethod = 'fit';


    %% ===== LEFT PANEL (Inputs) =====
    panel = uipanel(fig,...
        'Title','Blade Input Parameters',...
        'FontSize',14,...
        'Position',[20 20 300 560],...
        'BackgroundColor',[0 0 0],...
        'ForegroundColor','white');

    y = 500;
    spacing = 55;

    % Helper function for labels
    function createLabel(text,posY)
        uilabel(panel,...
            'Text',text,...
            'Position',[20 posY 200 22],...
            'FontColor','white',...
            'FontSize',12);
    end

    % Helper for numeric fields
    function field = createField(default,posY)
        field = uieditfield(panel,'numeric',...
            'Position',[190 posY 100 25],...
            'Value',default);
    end

    createLabel('Number of Mach Lines',y);
    nlField = createField(100,y); y=y-spacing;

    createLabel('Inlet Blade Angle (deg)',y);
    BetaField = createField(60,y); y=y-spacing;

    createLabel('Inlet Mach (M_i)',y);
    MiField = createField(3,y); y=y-spacing;

    createLabel('Pressure Surface Mach (M_L)',y);
    MLField = createField(2,y); y=y-spacing;

    createLabel('Suction Surface Mach (M_U)',y);
    MUField = createField(4,y); y=y-spacing;

    createLabel('Temperature (K)',y);
    TField = createField(300,y); y=y-spacing;

    createLabel('Gas Constant (R)',y);
    RField = createField(287,y); y=y-spacing;

    createLabel('Specific Heat Ratio (Y)',y);
    YField = createField(1.4,y); y=y-spacing;

    % Generate Button
    uibutton(panel,...
        'Text','Generate Blade',...
        'Position',[60 40 180 45],...
        'FontSize',14,...
        'BackgroundColor',[0.8 0.2 0.2],...
        'FontColor','white',...
        'ButtonPushedFcn',@(btn,event) generateBlade());

    function generateBlade()

        nl      = nlField.Value;
        Beta_i  = BetaField.Value;
        M_i     = MiField.Value;
        M_L     = MLField.Value;
        M_U     = MUField.Value;
        T       = TField.Value;
        R       = RField.Value;
        Y       = YField.Value;

        TREX_AUTO(nl, Beta_i, M_i, M_L, M_U, T, R, Y);
    end

    % -------------------------
    % Add bulleted points below the TREX image
    % Paste right after: img.ImageSource = 'TREX.png'; img.ScaleMethod = 'fit';
    % -------------------------
    panelPos = previewPanel.Position;    % [x y w h]
    pw = panelPos(3);
    ph = panelPos(4);
    
    % multi-line bulleted text (use \n for newlines)
    pointsText = sprintf(['\n• Specific Heat Ratio Y > 1.\n\n' ...
                          '• Upper(Suction) Surface Mach Number > Inlet Mach Number > Lower(Pressure) Surface Mach Number.\n\n' ...
                          '• M_U > M_i > M_L.\n\n' ...
                          '• Number of Mach lines > 2.\n\n' ...
                          '• Add more lines for a smoother curve (costs computation).\n\n' ...
                          '• CHECK YOUR BLADE PROFILES IN THE PROFILES FOLDER.\n\n']);
    
    % create a multi-line label centered horizontally near the bottom of previewPanel
    lblPoints = uilabel(previewPanel, ...
        'Text', pointsText, ...
        'FontColor',[1 1 1], ...
        'FontSize', 12, ...
        'HorizontalAlignment','left', ...
        'BackgroundColor','none', ...
        'WordWrap','on');
    
    % helper to (re)position the label in pixels according to panel size
    function placePointsLabel()
        panelPos = previewPanel.Position;
        pw = panelPos(3);
        ph = panelPos(4);
        w = round(1*pw);           % 85% of panel width
        h = round(0.5*ph);           % allocate ~22% of panel height for text block
        x = round((pw - w)/2);        % center horizontally
        y = round(0.08*ph);           % place near bottom (8% up from bottom)
        lblPoints.Position = [x, y, w, h];
        uistack(lblPoints,'top');     % ensure label is above the image
    end
    
    % call once now and whenever the previewPanel size changes
    placePointsLabel();

end