
package org.forester.archaeopteryx;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.List;

import javax.swing.ButtonGroup;
import javax.swing.JApplet;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

// Use like this:
// <applet archive="forester.jar"
// code="org.forester.archaeopteryx.ArchaeopteryxE.class"
// codebase="http://www.myserver.org/path/to/forester"
// width="600"
// height="500"
// alt="ArchaeopteryxE is not working on your system (requires at least Sun Java 1.5)!">
// <param name="url_of_tree_to_load"
// value="http://www.myserver.org/examples/data/apaf.xml">
// <param name="config_file"
// value="http://www.myserver.org/examples/config/config_file.txt">
// </applet>
public class ArchaeopteryxE extends JApplet implements ActionListener {

    private final static String  NAME             = "ArchaeopteryxE";
    private static final long    serialVersionUID = -1220055577935759443L;
    private Configuration        _configuration;
    private MainPanelApplets     _main_panel;
    private JMenuBar             _jmenubar;
    private JMenu                _options_jmenu;
    private JMenu                _font_size_menu;
    private JMenuItem            _super_tiny_fonts_mi;
    private JMenuItem            _tiny_fonts_mi;
    private JMenuItem            _small_fonts_mi;
    private JMenuItem            _medium_fonts_mi;
    private JMenuItem            _large_fonts_mi;
    private TextFrame            _textframe;
    private JMenu                _tools_menu;
    private JMenuItem            _taxcolor_item;
    private JMenuItem            _confcolor_item;
    private JMenuItem            _midpoint_root_item;
    private JMenu                _view_jmenu;
    private JMenuItem            _view_as_XML_item;
    private JMenuItem            _view_as_NH_item;
    private JMenuItem            _view_as_NHX_item;
    private JMenuItem            _view_as_nexus_item;
    private JMenu                _type_menu;
    private JCheckBoxMenuItem    _rectangular_type_cbmi;
    private JCheckBoxMenuItem    _triangular_type_cbmi;
    private JCheckBoxMenuItem    _curved_type_cbmi;
    private JCheckBoxMenuItem    _convex_type_cbmi;
    private JCheckBoxMenuItem    _euro_type_cbmi;
    private JCheckBoxMenuItem    _rounded_type_cbmi;
    private JCheckBoxMenuItem    _unrooted_type_cbmi;
    private JCheckBoxMenuItem    _circular_type_cbmi;
    private JMenuItem            _help_item;
    private JMenuItem            _about_item;
    private JMenu                _help_jmenu;
    private JMenuItem            _website_item;
    private JMenuItem            _phyloxml_website_item;
    private JMenuItem            _phyloxml_ref_item;
    private JMenuItem            _aptx_ref_item;
    private JMenuItem            _remove_branch_color_item;
    private JMenuItem            _infer_common_sn_names_item;
    private JCheckBoxMenuItem    _screen_antialias_cbmi;
    private JCheckBoxMenuItem    _background_gradient_cbmi;
    private JRadioButtonMenuItem _non_lined_up_cladograms_rbmi;
    private JRadioButtonMenuItem _uniform_cladograms_rbmi;
    private JRadioButtonMenuItem _ext_node_dependent_cladogram_rbmi;
    private Options              _options;
    private JMenuItem            _choose_font_mi;
    private JMenuItem            _switch_colors_mi;
    JCheckBoxMenuItem            _label_direction_cbmi;
    private JCheckBoxMenuItem    _show_node_boxes_cbmi;
    private JCheckBoxMenuItem    _show_scale_cbmi;
    private JCheckBoxMenuItem    _search_case_senstive_cbmi;
    private JCheckBoxMenuItem    _search_whole_words_only_cbmi;
    private JCheckBoxMenuItem    _inverse_search_result_cbmi;
    private JCheckBoxMenuItem    _show_overview_cbmi;
    private JMenuItem            _choose_minimal_confidence_mi;
    private JCheckBoxMenuItem    _show_branch_length_values_cbmi;
    private JMenuItem            _collapse_species_specific_subtrees;
    private JMenuItem            _overview_placment_mi;
    private JMenuItem            _collapse_below_threshold;
    private ButtonGroup          _radio_group_1;

    public void actionPerformed( final ActionEvent e ) {
        final Object o = e.getSource();
        if ( o == _midpoint_root_item ) {
            getMainPanel().getCurrentTreePanel().midpointRoot();
        }
        else if ( o == _taxcolor_item ) {
            getMainPanel().getCurrentTreePanel().taxColor();
        }
        else if ( o == _confcolor_item ) {
            getMainPanel().getCurrentTreePanel().confColor();
        }
        else if ( o == _infer_common_sn_names_item ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().inferCommonPartOfScientificNames();
            }
        }
        else if ( o == _collapse_species_specific_subtrees ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().collapseSpeciesSpecificSubtrees();
            }
        }
        else if ( o == _remove_branch_color_item ) {
            removeBranchColors();
        }
        else if ( o == _switch_colors_mi ) {
            switchColors();
        }
        else if ( o == _view_as_NH_item ) {
            viewAsNH();
        }
        else if ( o == _view_as_NHX_item ) {
            viewAsNHX();
        }
        else if ( o == _view_as_XML_item ) {
            viewAsXML();
        }
        else if ( o == _view_as_nexus_item ) {
            viewAsNexus();
        }
        else if ( o == _super_tiny_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setSuperTinyFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _tiny_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setTinyFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _small_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setSmallFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _medium_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setMediumFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _large_fonts_mi ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setLargeFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _choose_font_mi ) {
            chooseFont();
        }
        else if ( o == _choose_minimal_confidence_mi ) {
            chooseMinimalConfidence();
        }
        else if ( o == _overview_placment_mi ) {
            MainFrame.cycleOverview( getOptions(), getCurrentTreePanel() );
        }
        else if ( o == _show_node_boxes_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _non_lined_up_cladograms_rbmi ) {
            updateOptions( getOptions() );
            _main_panel.getControlPanel().showWhole();
        }
        else if ( o == _uniform_cladograms_rbmi ) {
            updateOptions( getOptions() );
            _main_panel.getControlPanel().showWhole();
        }
        else if ( o == _ext_node_dependent_cladogram_rbmi ) {
            updateOptions( getOptions() );
            _main_panel.getControlPanel().showWhole();
        }
        else if ( o == _search_case_senstive_cbmi ) {
            updateOptions( getOptions() );
            getMainPanel().getControlPanel().search();
        }
        else if ( o == _search_whole_words_only_cbmi ) {
            updateOptions( getOptions() );
            getMainPanel().getControlPanel().search();
        }
        else if ( o == _inverse_search_result_cbmi ) {
            updateOptions( getOptions() );
            getMainPanel().getControlPanel().search();
        }
        else if ( o == _show_scale_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _show_branch_length_values_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _label_direction_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _show_overview_cbmi ) {
            updateOptions( getOptions() );
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().updateOvSizes();
            }
        }
        else if ( ( o == _rectangular_type_cbmi ) || ( o == _triangular_type_cbmi ) || ( o == _curved_type_cbmi )
                || ( o == _convex_type_cbmi ) || ( o == _rounded_type_cbmi ) || ( o == _euro_type_cbmi )
                || ( o == _unrooted_type_cbmi ) || ( o == _circular_type_cbmi ) ) {
            typeChanged( o );
        }
        else if ( o == _screen_antialias_cbmi ) {
            updateOptions( getOptions() );
            setupScreenTextAntialias( getMainPanel().getTreePanels(), isScreenAntialias() );
        }
        else if ( o == _background_gradient_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _about_item ) {
            MainFrame.about();
        }
        else if ( o == _help_item ) {
            MainFrame.help( getConfiguration().getWebLinks() );
        }
        else if ( o == _website_item ) {
            try {
                Util.openWebsite( Constants.APTX_WEB_SITE, true, this );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _phyloxml_website_item ) {
            try {
                Util.openWebsite( Constants.PHYLOXML_WEB_SITE, true, this );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _aptx_ref_item ) {
            try {
                Util.openWebsite( Constants.APTX_REFERENCE_URL, true, this );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _phyloxml_ref_item ) {
            try {
                Util.openWebsite( Constants.PHYLOXML_REFERENCE_URL, true, this );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        repaint();
    }

    void buildFontSizeMenu() {
        _font_size_menu = MainFrame.createMenu( MainFrame.FONT_SIZE_MENU_LABEL, getConfiguration() );
        _font_size_menu.add( _super_tiny_fonts_mi = new JMenuItem( "Super tiny fonts" ) );
        _font_size_menu.add( _tiny_fonts_mi = new JMenuItem( "Tiny fonts" ) );
        _font_size_menu.add( _small_fonts_mi = new JMenuItem( "Small fonts" ) );
        _font_size_menu.add( _medium_fonts_mi = new JMenuItem( "Medium fonts" ) );
        _font_size_menu.add( _large_fonts_mi = new JMenuItem( "Large fonts" ) );
        customizeJMenuItem( _super_tiny_fonts_mi );
        customizeJMenuItem( _tiny_fonts_mi );
        customizeJMenuItem( _small_fonts_mi );
        customizeJMenuItem( _medium_fonts_mi );
        customizeJMenuItem( _large_fonts_mi );
        _jmenubar.add( _font_size_menu );
    }

    void buildHelpMenu() {
        _help_jmenu = MainFrame.createMenu( "Help", getConfiguration() );
        _help_jmenu.add( _help_item = new JMenuItem( "Help" ) );
        _help_jmenu.add( _website_item = new JMenuItem( "Archaeopteryx Home" ) );
        _aptx_ref_item = new JMenuItem( "Archaeopteryx Reference" );
        _help_jmenu.add( _phyloxml_website_item = new JMenuItem( "phyloXML Home" ) );
        _help_jmenu.add( _phyloxml_ref_item = new JMenuItem( "phyloXML Reference" ) );
        _help_jmenu.addSeparator();
        _help_jmenu.add( _about_item = new JMenuItem( "About" ) );
        customizeJMenuItem( _help_item );
        customizeJMenuItem( _website_item );
        customizeJMenuItem( _phyloxml_website_item );
        customizeJMenuItem( _aptx_ref_item );
        customizeJMenuItem( _phyloxml_ref_item );
        customizeJMenuItem( _about_item );
        _phyloxml_ref_item.setToolTipText( MainFrame.PHYLOXML_REF_TOOL_TIP );
        _aptx_ref_item.setToolTipText( MainFrame.APTX_REF_TOOL_TIP );
        _jmenubar.add( _help_jmenu );
    }

    void buildOptionsMenu() {
        _options_jmenu = MainFrame.createMenu( MainFrame.OPTIONS_HEADER, getConfiguration() );
        _options_jmenu.addChangeListener( new ChangeListener() {

            public void stateChanged( final ChangeEvent e ) {
                MainFrame.setOvPlacementColorChooseMenuItem( _overview_placment_mi, getCurrentTreePanel() );
                MainFrame.setTextColorChooseMenuItem( _switch_colors_mi, getCurrentTreePanel() );
                MainFrame
                        .setTextMinSupportMenuItem( _choose_minimal_confidence_mi, getOptions(), getCurrentTreePanel() );
                MainFrame.setTextForFontChooserMenuItem( _choose_font_mi, MainFrame
                        .createCurrentFontDesc( getMainPanel().getTreeFontSet() ) );
                MainFrame.updateOptionsMenuDependingOnPhylogenyType( getMainPanel(),
                                                                     _show_scale_cbmi,
                                                                     _show_branch_length_values_cbmi,
                                                                     _non_lined_up_cladograms_rbmi,
                                                                     _uniform_cladograms_rbmi,
                                                                     _ext_node_dependent_cladogram_rbmi,
                                                                     _label_direction_cbmi );
            }
        } );
        _options_jmenu.add( MainFrame.customizeMenuItemAsLabel( new JMenuItem( MainFrame.DISPLAY_SUBHEADER ),
                                                                getConfiguration() ) );
        _options_jmenu
                .add( _ext_node_dependent_cladogram_rbmi = new JRadioButtonMenuItem( MainFrame.NONUNIFORM_CLADOGRAMS_LABEL ) );
        _options_jmenu.add( _uniform_cladograms_rbmi = new JRadioButtonMenuItem( MainFrame.UNIFORM_CLADOGRAMS_LABEL ) );
        _options_jmenu
                .add( _non_lined_up_cladograms_rbmi = new JRadioButtonMenuItem( MainFrame.NON_LINED_UP_CLADOGRAMS_LABEL ) );
        _radio_group_1 = new ButtonGroup();
        _radio_group_1.add( _ext_node_dependent_cladogram_rbmi );
        _radio_group_1.add( _uniform_cladograms_rbmi );
        _radio_group_1.add( _non_lined_up_cladograms_rbmi );
        _options_jmenu.add( _show_node_boxes_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_NODE_BOXES_LABEL ) );
        _options_jmenu.add( _show_scale_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_SCALE_LABEL ) );
        _options_jmenu
                .add( _show_branch_length_values_cbmi = new JCheckBoxMenuItem( MainFrame.DISPLAY_BRANCH_LENGTH_VALUES_LABEL ) );
        _options_jmenu.add( _show_overview_cbmi = new JCheckBoxMenuItem( MainFrame.SHOW_OVERVIEW_LABEL ) );
        _options_jmenu.add( _label_direction_cbmi = new JCheckBoxMenuItem( MainFrame.LABEL_DIRECTION_LABEL ) );
        _label_direction_cbmi.setToolTipText( MainFrame.LABEL_DIRECTION_TIP );
        _options_jmenu.add( _screen_antialias_cbmi = new JCheckBoxMenuItem( MainFrame.SCREEN_ANTIALIAS_LABEL ) );
        _options_jmenu.add( _background_gradient_cbmi = new JCheckBoxMenuItem( MainFrame.BG_GRAD_LABEL ) );
        _options_jmenu.add( _choose_minimal_confidence_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _overview_placment_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _switch_colors_mi = new JMenuItem( "" ) );
        _options_jmenu.add( _choose_font_mi = new JMenuItem( "" ) );
        _options_jmenu.addSeparator();
        _options_jmenu.add( MainFrame.customizeMenuItemAsLabel( new JMenuItem( MainFrame.SEARCH_SUBHEADER ),
                                                                getConfiguration() ) );
        _options_jmenu
                .add( _search_case_senstive_cbmi = new JCheckBoxMenuItem( MainFrame.SEARCH_CASE_SENSITIVE_LABEL ) );
        _options_jmenu.add( _search_whole_words_only_cbmi = new JCheckBoxMenuItem( MainFrame.SEARCH_TERMS_ONLY_LABEL ) );
        _options_jmenu
                .add( _inverse_search_result_cbmi = new JCheckBoxMenuItem( MainFrame.INVERSE_SEARCH_RESULT_LABEL ) );
        customizeJMenuItem( _choose_font_mi );
        customizeJMenuItem( _choose_minimal_confidence_mi );
        customizeJMenuItem( _switch_colors_mi );
        customizeJMenuItem( _overview_placment_mi );
        customizeCheckBoxMenuItem( _show_node_boxes_cbmi, getOptions().isShowNodeBoxes() );
        customizeCheckBoxMenuItem( _label_direction_cbmi,
                                   getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL );
        customizeCheckBoxMenuItem( _screen_antialias_cbmi, getOptions().isAntialiasScreen() );
        customizeCheckBoxMenuItem( _background_gradient_cbmi, getOptions().isBackgroundColorGradient() );
        customizeCheckBoxMenuItem( _search_case_senstive_cbmi, getOptions().isSearchCaseSensitive() );
        customizeCheckBoxMenuItem( _show_scale_cbmi, getOptions().isShowScale() );
        customizeRadioButtonMenuItem( _non_lined_up_cladograms_rbmi,
                                      getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP );
        customizeRadioButtonMenuItem( _uniform_cladograms_rbmi,
                                      getOptions().getCladogramType() == CLADOGRAM_TYPE.TOTAL_NODE_SUM_DEP );
        customizeRadioButtonMenuItem( _ext_node_dependent_cladogram_rbmi,
                                      getOptions().getCladogramType() == CLADOGRAM_TYPE.EXT_NODE_SUM_DEP );
        customizeCheckBoxMenuItem( _show_branch_length_values_cbmi, getOptions().isShowBranchLengthValues() );
        customizeCheckBoxMenuItem( _show_overview_cbmi, getOptions().isShowOverview() );
        customizeCheckBoxMenuItem( _search_whole_words_only_cbmi, getOptions().isMatchWholeTermsOnly() );
        customizeCheckBoxMenuItem( _inverse_search_result_cbmi, getOptions().isInverseSearchResult() );
        _jmenubar.add( _options_jmenu );
    }

    void buildToolsMenu() {
        _tools_menu = MainFrame.createMenu( "Tools", getConfiguration() );
        _tools_menu.add( _confcolor_item = new JMenuItem( "Colorize Branches Depending on Confidence" ) );
        customizeJMenuItem( _confcolor_item );
        _tools_menu.add( _taxcolor_item = new JMenuItem( "Taxonomy Colorize Branches" ) );
        customizeJMenuItem( _taxcolor_item );
        _tools_menu.add( _remove_branch_color_item = new JMenuItem( "Delete Branch Colors" ) );
        _remove_branch_color_item.setToolTipText( "To delete branch color values from the current phylogeny." );
        customizeJMenuItem( _remove_branch_color_item );
        _tools_menu.addSeparator();
        _tools_menu.add( _midpoint_root_item = new JMenuItem( "Midpoint-Root" ) );
        customizeJMenuItem( _midpoint_root_item );
        _tools_menu.addSeparator();
        _tools_menu
                .add( _infer_common_sn_names_item = new JMenuItem( "Infer Common Parts of Internal Scientific Names" ) );
        customizeJMenuItem( _infer_common_sn_names_item );
        _tools_menu.add( _collapse_species_specific_subtrees = new JMenuItem( "Collapse Species-Specific Subtrees" ) );
        customizeJMenuItem( _collapse_species_specific_subtrees );
        if ( !ForesterConstants.RELEASE ) {
            _tools_menu
                    .add( _collapse_below_threshold = new JMenuItem( "Collapse Branches with Confidence Below Threshold" ) );
            customizeJMenuItem( _collapse_below_threshold );
        }
        _jmenubar.add( _tools_menu );
    }

    void buildTypeMenu() {
        _type_menu = MainFrame.createMenu( MainFrame.TYPE_MENU_HEADER, getConfiguration() );
        _type_menu.add( _rectangular_type_cbmi = new JCheckBoxMenuItem( MainFrame.RECTANGULAR_TYPE_CBMI_LABEL ) );
        _type_menu.add( _euro_type_cbmi = new JCheckBoxMenuItem( MainFrame.EURO_TYPE_CBMI_LABEL ) );
        _type_menu.add( _rounded_type_cbmi = new JCheckBoxMenuItem( MainFrame.ROUNDED_TYPE_CBMI_LABEL ) );
        _type_menu.add( _curved_type_cbmi = new JCheckBoxMenuItem( MainFrame.CURVED_TYPE_CBMI_LABEL ) );
        _type_menu.add( _triangular_type_cbmi = new JCheckBoxMenuItem( MainFrame.TRIANGULAR_TYPE_CBMI_LABEL ) );
        _type_menu.add( _convex_type_cbmi = new JCheckBoxMenuItem( MainFrame.CONVEX_TYPE_CBMI_LABEL ) );
        _type_menu.add( _unrooted_type_cbmi = new JCheckBoxMenuItem( MainFrame.UNROOTED_TYPE_CBMI_LABEL ) );
        _type_menu.add( _circular_type_cbmi = new JCheckBoxMenuItem( MainFrame.CIRCULAR_TYPE_CBMI_LABEL ) );
        customizeCheckBoxMenuItem( _rectangular_type_cbmi, false );
        customizeCheckBoxMenuItem( _triangular_type_cbmi, false );
        customizeCheckBoxMenuItem( _euro_type_cbmi, false );
        customizeCheckBoxMenuItem( _rounded_type_cbmi, false );
        customizeCheckBoxMenuItem( _curved_type_cbmi, false );
        customizeCheckBoxMenuItem( _convex_type_cbmi, false );
        customizeCheckBoxMenuItem( _unrooted_type_cbmi, false );
        customizeCheckBoxMenuItem( _circular_type_cbmi, false );
        _unrooted_type_cbmi.setToolTipText( MainFrame.USE_MOUSEWHEEL_SHIFT_TO_ROTATE );
        _circular_type_cbmi.setToolTipText( MainFrame.USE_MOUSEWHEEL_SHIFT_TO_ROTATE );
        initializeTypeMenu( getOptions() );
        _jmenubar.add( _type_menu );
    }

    void buildViewMenu() {
        _view_jmenu = MainFrame.createMenu( "View as Text", getConfiguration() );
        _view_jmenu.add( _view_as_XML_item = new JMenuItem( "View as phyloXML" ) );
        _view_jmenu.add( _view_as_NH_item = new JMenuItem( "View as Newick" ) );
        _view_jmenu.add( _view_as_NHX_item = new JMenuItem( "View as NHX" ) );
        _view_jmenu.add( _view_as_nexus_item = new JMenuItem( "View as Nexus" ) );
        customizeJMenuItem( _view_as_NH_item );
        customizeJMenuItem( _view_as_NHX_item );
        customizeJMenuItem( _view_as_XML_item );
        customizeJMenuItem( _view_as_nexus_item );
        _jmenubar.add( _view_jmenu );
    }

    private void chooseFont() {
        final FontChooser fc = new FontChooser();
        fc.setFont( getMainPanel().getTreeFontSet().getLargeFont() );
        fc.showDialog( this, "Select the Base Font" );
        getMainPanel().getTreeFontSet().setBaseFont( fc.getFont() );
    }

    private void chooseMinimalConfidence() {
        final String s = ( String ) JOptionPane
                .showInputDialog( this,
                                  "Please the minimum for confidence values to be displayed.\n" + "[current value: "
                                          + getOptions().getMinConfidenceValue() + "]\n",
                                  "Minimal Confidence Value",
                                  JOptionPane.QUESTION_MESSAGE,
                                  null,
                                  null,
                                  getOptions().getMinConfidenceValue() );
        if ( !ForesterUtil.isEmpty( s ) ) {
            boolean success = true;
            double m = 0.0;
            final String m_str = s.trim();
            if ( !ForesterUtil.isEmpty( m_str ) ) {
                try {
                    m = Double.parseDouble( m_str );
                }
                catch ( final Exception ex ) {
                    success = false;
                }
            }
            else {
                success = false;
            }
            if ( success && ( m >= 0.0 ) ) {
                getOptions().setMinConfidenceValue( m );
            }
        }
    }

    void customizeCheckBoxMenuItem( final JCheckBoxMenuItem item, final boolean is_selected ) {
        if ( item != null ) {
            item.setFont( MainFrame.menu_font );
            if ( !getConfiguration().isUseNativeUI() ) {
                item.setBackground( Constants.MENU_BACKGROUND_COLOR_DEFAULT );
                item.setForeground( Constants.MENU_TEXT_COLOR_DEFAULT );
            }
            item.setSelected( is_selected );
            item.addActionListener( this );
        }
    }

    void customizeJMenuItem( final JMenuItem jmi ) {
        jmi.setFont( MainFrame.menu_font );
        if ( !getConfiguration().isUseNativeUI() ) {
            jmi.setBackground( Constants.MENU_BACKGROUND_COLOR_DEFAULT );
            jmi.setForeground( Constants.MENU_TEXT_COLOR_DEFAULT );
        }
        jmi.addActionListener( this );
    }

    private void customizeRadioButtonMenuItem( final JRadioButtonMenuItem item, final boolean is_selected ) {
        if ( item != null ) {
            item.setFont( MainFrame.menu_font );
            if ( !getConfiguration().isUseNativeUI() ) {
                item.setBackground( Constants.MENU_BACKGROUND_COLOR_DEFAULT );
                item.setForeground( Constants.MENU_TEXT_COLOR_DEFAULT );
            }
            item.setSelected( is_selected );
            item.addActionListener( this );
        }
    }

    @Override
    public void destroy() {
        Util.printAppletMessage( NAME, "going to be destroyed " );
        removeTextFrame();
        if ( getMainPanel() != null ) {
            getMainPanel().terminate();
        }
    }

    Configuration getConfiguration() {
        return _configuration;
    }

    TreePanel getCurrentTreePanel() {
        return getMainPanel().getCurrentTreePanel();
    }

    JCheckBoxMenuItem getlabelDirectionCbmi() {
        return _label_direction_cbmi;
    }

    private MainPanel getMainPanel() {
        return _main_panel;
    }

    public Options getOptions() {
        return _options;
    }

    Options getOtions() {
        return _options;
    }

    @Override
    public void init() {
        final String config_filename = getParameter( Constants.APPLET_PARAM_NAME_FOR_CONFIG_FILE_URL );
        Util.printAppletMessage( NAME, "URL for configuration file is: " + config_filename );
        final Configuration configuration = new Configuration( config_filename, true, true );
        setConfiguration( configuration );
        setOptions( Options.createInstance( configuration ) );
        setupUI();
        URL phys_url = null;
        Phylogeny[] phys = null;
        final String phys_url_string = getParameter( Constants.APPLET_PARAM_NAME_FOR_URL_OF_TREE_TO_LOAD );
        Util.printAppletMessage( NAME, "URL for phylogenies is " + phys_url_string );
        // Get URL to tree file
        if ( phys_url_string != null ) {
            try {
                phys_url = new URL( phys_url_string );
            }
            catch ( final Exception e ) {
                ForesterUtil.printErrorMessage( NAME, "error: " + e );
                e.printStackTrace();
                JOptionPane.showMessageDialog( this, NAME + ": Could not create URL from: \"" + phys_url_string
                        + "\"\nException: " + e, "Failed to create URL", JOptionPane.ERROR_MESSAGE );
            }
        }
        // Load the tree from URL
        if ( phys_url != null ) {
            try {
                phys = Util.readPhylogeniesFromUrl( phys_url, getConfiguration().isValidatePhyloXmlAgainstSchema() );
            }
            catch ( final Exception e ) {
                ForesterUtil.printErrorMessage( NAME, e.toString() );
                e.printStackTrace();
                JOptionPane.showMessageDialog( this,
                                               NAME + ": Failed to read phylogenies: " + "\nException: " + e,
                                               "Failed to read phylogenies",
                                               JOptionPane.ERROR_MESSAGE );
            }
        }
        if ( ( phys == null ) || ( phys.length < 1 ) ) {
            ForesterUtil.printErrorMessage( NAME, "phylogenies from [" + phys_url + "] are null or empty" );
            JOptionPane.showMessageDialog( this,
                                           NAME + ": phylogenies from [" + phys_url + "] are null or empty",
                                           "Failed to read phylogenies",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        else {
            Util.printAppletMessage( NAME, "loaded " + phys.length + " phylogenies from: " + phys_url );
        }
        setVisible( false );
        setMainPanel( new MainPanelApplets( getConfiguration(), this ) );
        _jmenubar = new JMenuBar();
        if ( !getConfiguration().isHideControlPanelAndMenubar() ) {
            if ( !getConfiguration().isUseNativeUI() ) {
                _jmenubar.setBackground( Constants.MENU_BACKGROUND_COLOR_DEFAULT );
            }
            buildToolsMenu();
            buildViewMenu();
            buildFontSizeMenu();
            buildOptionsMenu();
            buildTypeMenu();
            buildHelpMenu();
            setJMenuBar( _jmenubar );
        }
        final Container contentpane = getContentPane();
        contentpane.setLayout( new BorderLayout() );
        contentpane.add( getMainPanel(), BorderLayout.CENTER );
        addComponentListener( new ComponentAdapter() {

            @Override
            public void componentResized( final ComponentEvent e ) {
                if ( getMainPanel().getCurrentTreePanel() != null ) {
                    getMainPanel().getCurrentTreePanel().setParametersForPainting( getMainPanel().getCurrentTreePanel()
                                                                                           .getWidth(),
                                                                                   getMainPanel().getCurrentTreePanel()
                                                                                           .getHeight(),
                                                                                   false );
                }
            }
        } );
        if ( getConfiguration().isUseTabbedDisplay() ) {
            Util.printAppletMessage( NAME, "using tabbed display" );
            Util.addPhylogeniesToTabs( phys,
                                       new File( phys_url.getFile() ).getName(),
                                       phys_url.toString(),
                                       getConfiguration(),
                                       getMainPanel() );
        }
        else {
            Util.printAppletMessage( NAME, "not using tabbed display" );
            Util.addPhylogenyToPanel( phys, getConfiguration(), getMainPanel() );
        }
        validate();
        setName( NAME );
        getMainPanel().getControlPanel().showWholeAll();
        getMainPanel().getControlPanel().showWhole();
        System.gc();
        Util.printAppletMessage( NAME, "successfully initialized" );
        setVisible( true );
    }

    void initializeTypeMenu( final Options options ) {
        setTypeMenuToAllUnselected();
        try {
            switch ( options.getPhylogenyGraphicsType() ) {
                case CONVEX:
                    _convex_type_cbmi.setSelected( true );
                    break;
                case CURVED:
                    _curved_type_cbmi.setSelected( true );
                    break;
                case EURO_STYLE:
                    _euro_type_cbmi.setSelected( true );
                    break;
                case ROUNDED:
                    _rounded_type_cbmi.setSelected( true );
                    break;
                case TRIANGULAR:
                    _triangular_type_cbmi.setSelected( true );
                    break;
                case UNROOTED:
                    _unrooted_type_cbmi.setSelected( true );
                    break;
                case CIRCULAR:
                    _circular_type_cbmi.setSelected( true );
                    break;
                default:
                    _rectangular_type_cbmi.setSelected( true );
                    break;
            }
        }
        catch ( final NullPointerException np ) {
            // In all likelihood, this is caused by menu-less display.
        }
    }

    private boolean isScreenAntialias() {
        return true;
    }

    private void removeBranchColors() {
        if ( getMainPanel().getCurrentPhylogeny() != null ) {
            Util.removeBranchColors( getMainPanel().getCurrentPhylogeny() );
        }
    }

    void removeTextFrame() {
        if ( _textframe != null ) {
            _textframe.close();
            _textframe = null;
        }
    }

    void setConfiguration( final Configuration configuration ) {
        _configuration = configuration;
    }

    private void setMainPanel( final MainPanelApplets main_panel ) {
        _main_panel = main_panel;
    }

    void setOptions( final Options options ) {
        _options = options;
    }

    void setSelectedTypeInTypeMenu( final PHYLOGENY_GRAPHICS_TYPE type ) {
        setTypeMenuToAllUnselected();
        try {
            switch ( type ) {
                case CIRCULAR:
                    _circular_type_cbmi.setSelected( true );
                    break;
                case CONVEX:
                    _convex_type_cbmi.setSelected( true );
                    break;
                case CURVED:
                    _curved_type_cbmi.setSelected( true );
                    break;
                case EURO_STYLE:
                    _euro_type_cbmi.setSelected( true );
                    break;
                case ROUNDED:
                    _rounded_type_cbmi.setSelected( true );
                    break;
                case RECTANGULAR:
                    _rectangular_type_cbmi.setSelected( true );
                    break;
                case TRIANGULAR:
                    _triangular_type_cbmi.setSelected( true );
                    break;
                case UNROOTED:
                    _unrooted_type_cbmi.setSelected( true );
                    break;
                default:
                    throw new IllegalArgumentException( "unknown type: " + type );
            }
        }
        catch ( final NullPointerException np ) {
            // In all likelihood, this is caused by menu-less display.
        }
    }

    void setTypeMenuToAllUnselected() {
        if ( _convex_type_cbmi != null ) {
            _convex_type_cbmi.setSelected( false );
        }
        if ( _curved_type_cbmi != null ) {
            _curved_type_cbmi.setSelected( false );
        }
        if ( _euro_type_cbmi != null ) {
            _euro_type_cbmi.setSelected( false );
        }
        if ( _rounded_type_cbmi != null ) {
            _rounded_type_cbmi.setSelected( false );
        }
        if ( _triangular_type_cbmi != null ) {
            _triangular_type_cbmi.setSelected( false );
        }
        if ( _rectangular_type_cbmi != null ) {
            _rectangular_type_cbmi.setSelected( false );
        }
        if ( _unrooted_type_cbmi != null ) {
            _unrooted_type_cbmi.setSelected( false );
        }
        if ( _circular_type_cbmi != null ) {
            _circular_type_cbmi.setSelected( false );
        }
    }

    private void setupUI() {
        try {
            if ( getConfiguration().isUseNativeUI() ) {
                UIManager.setLookAndFeel( UIManager.getSystemLookAndFeelClassName() );
            }
            else {
                UIManager.setLookAndFeel( UIManager.getCrossPlatformLookAndFeelClassName() );
            }
        }
        catch ( final UnsupportedLookAndFeelException e ) {
            Util.dieWithSystemError( "UnsupportedLookAndFeelException: " + e.toString() );
        }
        catch ( final ClassNotFoundException e ) {
            Util.dieWithSystemError( "ClassNotFoundException: " + e.toString() );
        }
        catch ( final InstantiationException e ) {
            Util.dieWithSystemError( "InstantiationException: " + e.toString() );
        }
        catch ( final IllegalAccessException e ) {
            Util.dieWithSystemError( "IllegalAccessException: " + e.toString() );
        }
        catch ( final Exception e ) {
            Util.dieWithSystemError( e.toString() );
        }
    }

    @Override
    public void start() {
        if ( getMainPanel() != null ) {
            getMainPanel().validate();
        }
        requestFocus();
        requestFocusInWindow();
        requestFocus();
        Util.printAppletMessage( NAME, "started" );
    }

    void switchColors() {
        final TreeColorSet colorset = getMainPanel().getCurrentTreePanel().getTreeColorSet();
        final ColorSchemeChooser csc = new ColorSchemeChooser( getMainPanel(), colorset );
        csc.setVisible( true );
        getMainPanel().setTreeColorSet( colorset );
    }

    void typeChanged( final Object o ) {
        updateTypeCheckboxes( getOptions(), o );
        updateOptions( getOptions() );
        if ( getCurrentTreePanel() != null ) {
            final PHYLOGENY_GRAPHICS_TYPE previous_type = getCurrentTreePanel().getPhylogenyGraphicsType();
            final PHYLOGENY_GRAPHICS_TYPE new_type = getOptions().getPhylogenyGraphicsType();
            if ( ( ( previous_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && ( new_type != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) )
                    || ( ( previous_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) && ( new_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) )
                    || ( ( previous_type != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && ( new_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) )
                    || ( ( previous_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) && ( new_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) ) {
                getCurrentTreePanel().getControlPanel().showWhole();
            }
            if ( getCurrentTreePanel().isPhyHasBranchLengths() && ( new_type != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                getCurrentTreePanel().getControlPanel().setDrawPhylogramEnabled( true );
            }
            else {
                getCurrentTreePanel().getControlPanel().setDrawPhylogramEnabled( false );
            }
            getCurrentTreePanel().setPhylogenyGraphicsType( getOptions().getPhylogenyGraphicsType() );
            MainFrame.updateScreenTextAntialias( getMainPanel().getTreePanels() );
        }
    }

    void updateOptions( final Options options ) {
        options.setAntialiasScreen( ( _screen_antialias_cbmi != null ) && _screen_antialias_cbmi.isSelected() );
        options.setBackgroundColorGradient( ( _background_gradient_cbmi != null )
                && _background_gradient_cbmi.isSelected() );
        options.setShowNodeBoxes( ( _show_node_boxes_cbmi != null ) && _show_node_boxes_cbmi.isSelected() );
        if ( ( _non_lined_up_cladograms_rbmi != null ) && ( _non_lined_up_cladograms_rbmi.isSelected() ) ) {
            options.setCladogramType( CLADOGRAM_TYPE.NON_LINED_UP );
        }
        else if ( ( _uniform_cladograms_rbmi != null ) && ( _uniform_cladograms_rbmi.isSelected() ) ) {
            options.setCladogramType( CLADOGRAM_TYPE.TOTAL_NODE_SUM_DEP );
        }
        else if ( ( _ext_node_dependent_cladogram_rbmi != null ) && ( _ext_node_dependent_cladogram_rbmi.isSelected() ) ) {
            options.setCladogramType( CLADOGRAM_TYPE.EXT_NODE_SUM_DEP );
        }
        options.setSearchCaseSensitive( ( _search_case_senstive_cbmi != null )
                && _search_case_senstive_cbmi.isSelected() );
        if ( ( _show_scale_cbmi != null ) && _show_scale_cbmi.isEnabled() ) {
            options.setShowScale( _show_scale_cbmi.isSelected() );
        }
        if ( _label_direction_cbmi != null ) {
            if ( _label_direction_cbmi.isSelected() ) {
                options.setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
            }
            else {
                options.setNodeLabelDirection( NODE_LABEL_DIRECTION.HORIZONTAL );
            }
        }
        options.setShowOverview( ( _show_overview_cbmi != null ) && _show_overview_cbmi.isSelected() );
        if ( ( _show_branch_length_values_cbmi != null ) && _show_branch_length_values_cbmi.isEnabled() ) {
            options.setShowBranchLengthValues( _show_branch_length_values_cbmi.isSelected() );
        }
        options.setMatchWholeTermsOnly( ( _search_whole_words_only_cbmi != null )
                && _search_whole_words_only_cbmi.isSelected() );
        options.setInverseSearchResult( ( _inverse_search_result_cbmi != null )
                && _inverse_search_result_cbmi.isSelected() );
        if ( ( _rectangular_type_cbmi != null ) && _rectangular_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        }
        else if ( ( _triangular_type_cbmi != null ) && _triangular_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR );
        }
        else if ( ( _curved_type_cbmi != null ) && _curved_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CURVED );
        }
        else if ( ( _convex_type_cbmi != null ) && _convex_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CONVEX );
        }
        else if ( ( _euro_type_cbmi != null ) && _euro_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE );
        }
        else if ( ( _rounded_type_cbmi != null ) && _rounded_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.ROUNDED );
        }
        else if ( ( _unrooted_type_cbmi != null ) && _unrooted_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
        }
        else if ( ( _circular_type_cbmi != null ) && _circular_type_cbmi.isSelected() ) {
            options.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        }
    }

    void updateTypeCheckboxes( final Options options, final Object o ) {
        setTypeMenuToAllUnselected();
        ( ( JCheckBoxMenuItem ) o ).setSelected( true );
    }

    void viewAsNexus() {
        removeTextFrame();
        if ( ( getMainPanel().getCurrentPhylogeny() == null ) || getMainPanel().getCurrentPhylogeny().isEmpty()
                || ( getMainPanel().getCurrentPhylogeny().getNumberOfExternalNodes() > 10000 ) ) {
            return;
        }
        _textframe = TextFrame.instantiate( getMainPanel().getCurrentPhylogeny().toNexus() );
    }

    void viewAsNH() {
        removeTextFrame();
        if ( ( getMainPanel().getCurrentPhylogeny() == null ) || getMainPanel().getCurrentPhylogeny().isEmpty()
                || ( getMainPanel().getCurrentPhylogeny().getNumberOfExternalNodes() > 10000 ) ) {
            return;
        }
        _textframe = TextFrame.instantiate( getMainPanel().getCurrentPhylogeny().toNewHampshire( false ) );
    }

    void viewAsNHX() {
        removeTextFrame();
        if ( ( getMainPanel().getCurrentPhylogeny() == null ) || getMainPanel().getCurrentPhylogeny().isEmpty()
                || ( getMainPanel().getCurrentPhylogeny().getNumberOfExternalNodes() > 10000 ) ) {
            return;
        }
        _textframe = TextFrame.instantiate( getMainPanel().getCurrentPhylogeny().toNewHampshireX() );
    }

    void viewAsXML() {
        removeTextFrame();
        if ( ( getMainPanel().getCurrentPhylogeny() == null ) || getMainPanel().getCurrentPhylogeny().isEmpty()
                || ( getMainPanel().getCurrentPhylogeny().getNumberOfExternalNodes() > 10000 ) ) {
            return;
        }
        _textframe = TextFrame.instantiate( getMainPanel().getCurrentPhylogeny().toPhyloXML( 0 ) );
    }

    static void setupScreenTextAntialias( final List<TreePanel> treepanels, final boolean antialias ) {
        for( final TreePanel tree_panel : treepanels ) {
            tree_panel.setTextAntialias();
        }
    }
}