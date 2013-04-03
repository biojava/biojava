// $Id: MainFrame.java,v 1.67 2009/12/09 20:13:50 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2010 Christian M. Zmasek
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek at yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx;

import java.awt.Container;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import javax.swing.JApplet;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.webservices.PhylogeniesWebserviceClient;
import org.forester.archaeopteryx.webservices.WebserviceUtil;
import org.forester.archaeopteryx.webservices.WebservicesManager;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.tol.TolParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public abstract class MainFrame extends JFrame implements ActionListener {

    static final String       USE_MOUSEWHEEL_SHIFT_TO_ROTATE     = "In this display type, use mousewheel + Shift to rotate [or A and S]";
    static final String       PHYLOXML_REF_TOOL_TIP              = Constants.PHYLOXML_REFERENCE;                                         //TODO //FIXME
    static final String       APTX_REF_TOOL_TIP                  = Constants.APTX_REFERENCE;
    private static final long serialVersionUID                   = 3655000897845508358L;
    final static Font         menu_font                          = new Font( Configuration.getDefaultFontFamilyName(),
                                                                             Font.PLAIN,
                                                                             10 );
    static final String       TYPE_MENU_HEADER                   = "Type";
    static final String       RECTANGULAR_TYPE_CBMI_LABEL        = "Rectangular";
    static final String       EURO_TYPE_CBMI_LABEL               = "Euro Type";
    static final String       CURVED_TYPE_CBMI_LABEL             = "Curved";
    static final String       TRIANGULAR_TYPE_CBMI_LABEL         = "Triangular";
    static final String       CONVEX_TYPE_CBMI_LABEL             = "Convex";
    static final String       ROUNDED_TYPE_CBMI_LABEL            = "Rounded";
    static final String       UNROOTED_TYPE_CBMI_LABEL           = "Unrooted (alpha)";                                                   //TODO
    static final String       CIRCULAR_TYPE_CBMI_LABEL           = "Circular (alpha)";                                                   //TODO
    static final String       OPTIONS_HEADER                     = "Options";
    static final String       SEARCH_SUBHEADER                   = "Search:";
    static final String       DISPLAY_SUBHEADER                  = "Display:";
    static final String       SEARCH_TERMS_ONLY_LABEL            = "Match Complete Terms Only";
    static final String       SEARCH_CASE_SENSITIVE_LABEL        = "Case Sensitive";
    static final String       INVERSE_SEARCH_RESULT_LABEL        = "Negate Result";
    static final String       DISPLAY_BRANCH_LENGTH_VALUES_LABEL = "Display Branch Length Values";
    static final String       DISPLAY_SCALE_LABEL                = "Display Scale";
    static final String       NON_LINED_UP_CLADOGRAMS_LABEL      = "Non-Lined Up Cladograms";
    static final String       UNIFORM_CLADOGRAMS_LABEL           = "Total Node Sum Dependent Cladograms";
    static final String       LABEL_DIRECTION_LABEL              = "Radial Labels";
    static final String       LABEL_DIRECTION_TIP                = "To use radial node labels in radial and unrooted display types";
    static final String       SCREEN_ANTIALIAS_LABEL             = "Antialias";
    static final String       BG_GRAD_LABEL                      = "Background Color Gradient";
    static final String       DISPLAY_NODE_BOXES_LABEL           = "Display Node Boxes";
    static final String       SHOW_OVERVIEW_LABEL                = "Show Overview";
    static final String       FONT_SIZE_MENU_LABEL               = "Font Size";
    static final String       NONUNIFORM_CLADOGRAMS_LABEL        = "External Node Sum Dependent Cladograms";
    JMenuBar                  _jmenubar;
    JMenu                     _file_jmenu;
    JMenu                     _tools_menu;
    JMenu                     _view_jmenu;
    JMenu                     _options_jmenu;
    JMenu                     _font_size_menu;
    JMenu                     _help_jmenu;
    JMenuItem[]               _load_phylogeny_from_webservice_menu_items;
    // file menu:
    JMenuItem                 _open_item;
    JMenuItem                 _open_url_item;
    JMenuItem                 _save_item;
    JMenuItem                 _save_all_item;
    JMenuItem                 _close_item;
    JMenuItem                 _exit_item;
    JMenuItem                 _new_item;
    // tools menu:
    JMenuItem                 _midpoint_root_item;
    JMenuItem                 _taxcolor_item;
    JMenuItem                 _confcolor_item;
    JMenuItem                 _infer_common_sn_names_item;
    JMenuItem                 _collapse_species_specific_subtrees;
    JMenuItem                 _collapse_below_threshold;                                                                                 //TODO implememt me
    // font size menu:
    JMenuItem                 _super_tiny_fonts_item;
    JMenuItem                 _tiny_fonts_item;
    JMenuItem                 _small_fonts_item;
    JMenuItem                 _medium_fonts_item;
    JMenuItem                 _large_fonts_item;
    // options menu:
    // _  screen and print
    JMenuItem                 _choose_font_mi;
    JMenuItem                 _switch_colors_mi;
    JCheckBoxMenuItem         _label_direction_cbmi;
    // _  screen display
    JCheckBoxMenuItem         _screen_antialias_cbmi;
    JCheckBoxMenuItem         _background_gradient_cbmi;
    JCheckBoxMenuItem         _show_node_boxes_cbmi;
    JRadioButtonMenuItem      _non_lined_up_cladograms_rbmi;
    JRadioButtonMenuItem      _uniform_cladograms_rbmi;
    JRadioButtonMenuItem      _ext_node_dependent_cladogram_rbmi;
    JCheckBoxMenuItem         _show_branch_length_values_cbmi;
    JCheckBoxMenuItem         _show_scale_cbmi;                                                                                          //TODO fix me
    JCheckBoxMenuItem         _show_overview_cbmi;
    JMenuItem                 _overview_placment_mi;
    JMenuItem                 _choose_minimal_confidence_mi;
    // _  print
    JCheckBoxMenuItem         _graphics_export_visible_only_cbmi;
    JCheckBoxMenuItem         _antialias_print_cbmi;
    JCheckBoxMenuItem         _print_black_and_white_cbmi;
    JCheckBoxMenuItem         _print_using_actual_size_cbmi;
    JCheckBoxMenuItem         _graphics_export_using_actual_size_cbmi;
    JMenuItem                 _print_size_mi;
    JMenuItem                 _choose_pdf_width_mi;
    // _  parsing
    JCheckBoxMenuItem         _internal_number_are_confidence_for_nh_parsing_cbmi;
    JCheckBoxMenuItem         _extract_pfam_style_tax_codes_cbmi;
    JCheckBoxMenuItem         _replace_underscores_cbmi;
    JCheckBoxMenuItem         _ignore_quotes_in_nh_parsing_cbmi;
    // _  search
    JCheckBoxMenuItem         _search_case_senstive_cbmi;
    JCheckBoxMenuItem         _search_whole_words_only_cbmi;
    JCheckBoxMenuItem         _inverse_search_result_cbmi;
    // type menu:
    JMenu                     _type_menu;
    JCheckBoxMenuItem         _rectangular_type_cbmi;
    JCheckBoxMenuItem         _triangular_type_cbmi;
    JCheckBoxMenuItem         _curved_type_cbmi;
    JCheckBoxMenuItem         _convex_type_cbmi;
    JCheckBoxMenuItem         _euro_type_cbmi;
    JCheckBoxMenuItem         _rounded_type_cbmi;
    JCheckBoxMenuItem         _unrooted_type_cbmi;
    JCheckBoxMenuItem         _circular_type_cbmi;
    // view as text menu:
    JMenuItem                 _view_as_NH_item;
    JMenuItem                 _view_as_NHX_item;
    JMenuItem                 _view_as_XML_item;
    JMenuItem                 _view_as_nexus_item;
    // help menu:
    JMenuItem                 _about_item;
    JMenuItem                 _help_item;
    JMenuItem                 _website_item;
    JMenuItem                 _phyloxml_website_item;
    JMenuItem                 _phyloxml_ref_item;
    JMenuItem                 _aptx_ref_item;
    // Handy pointers to child components:
    MainPanel                 _mainpanel;
    Container                 _contentpane;
    TextFrame                 _textframe;
    Configuration             _configuration;
    JMenuItem                 _remove_branch_color_item;
    Options                   _options;

    MainFrame() {
        // Empty constructor.
    }

    /**
     * Action performed.
     */
    public void actionPerformed( final ActionEvent e ) {
        final Object o = e.getSource();
        boolean is_applet = false;
        JApplet applet = null;
        if ( getCurrentTreePanel() != null ) {
            is_applet = getCurrentTreePanel().isApplet();
            if ( is_applet ) {
                applet = getCurrentTreePanel().obtainApplet();
            }
        }
        if ( o == _open_url_item ) {
            readPhylogeniesFromURL();
        }
        else if ( o == _exit_item ) {
            close();
        }
        else if ( o == _taxcolor_item ) {
            taxColor();
        }
        else if ( o == _confcolor_item ) {
            confColor();
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
        else if ( o == _midpoint_root_item ) {
            midpointRoot();
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
        else if ( o == _super_tiny_fonts_item ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setSuperTinyFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _tiny_fonts_item ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setTinyFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _small_fonts_item ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setSmallFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _medium_fonts_item ) {
            if ( getCurrentTreePanel() != null ) {
                getCurrentTreePanel().setMediumFonts();
                getCurrentTreePanel().repaint();
            }
        }
        else if ( o == _large_fonts_item ) {
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
        else if ( o == _screen_antialias_cbmi ) {
            updateOptions( getOptions() );
            updateScreenTextAntialias( getMainPanel().getTreePanels() );
        }
        else if ( o == _background_gradient_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _show_node_boxes_cbmi ) {
            updateOptions( getOptions() );
        }
        else if ( o == _non_lined_up_cladograms_rbmi ) {
            updateOptions( getOptions() );
            showWhole();
        }
        else if ( o == _uniform_cladograms_rbmi ) {
            updateOptions( getOptions() );
            showWhole();
        }
        else if ( o == _ext_node_dependent_cladogram_rbmi ) {
            updateOptions( getOptions() );
            showWhole();
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
                || ( o == _convex_type_cbmi ) || ( o == _euro_type_cbmi ) || ( o == _rounded_type_cbmi )
                || ( o == _unrooted_type_cbmi ) || ( o == _circular_type_cbmi ) ) {
            typeChanged( o );
        }
        else if ( o == _about_item ) {
            about();
        }
        else if ( o == _help_item ) {
            help( getConfiguration().getWebLinks() );
        }
        else if ( o == _website_item ) {
            try {
                Util.openWebsite( Constants.APTX_WEB_SITE, is_applet, applet );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _phyloxml_website_item ) {
            try {
                Util.openWebsite( Constants.PHYLOXML_WEB_SITE, is_applet, applet );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _aptx_ref_item ) {
            try {
                Util.openWebsite( Constants.APTX_REFERENCE_URL, is_applet, applet );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else if ( o == _phyloxml_ref_item ) {
            try {
                Util.openWebsite( Constants.PHYLOXML_REFERENCE_URL, is_applet, applet );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.printErrorMessage( Constants.PRG_NAME, e1.toString() );
            }
        }
        else {
            for( int i = 0; i < _load_phylogeny_from_webservice_menu_items.length; ++i ) {
                if ( o == _load_phylogeny_from_webservice_menu_items[ i ] ) {
                    readPhylogeniesFromWebservice( i );
                }
            }
        }
        _contentpane.repaint();
    }

    void activateSaveAllIfNeeded() {
        if ( ( getMainPanel().getTabbedPane() != null ) && ( getMainPanel().getTabbedPane().getTabCount() > 1 ) ) {
            _save_all_item.setEnabled( true );
        }
        else {
            _save_all_item.setEnabled( false );
        }
    }

    void buildFileMenu() {
        _file_jmenu = createMenu( "File", getConfiguration() );
        _file_jmenu.add( _open_url_item = new JMenuItem( "Read tree from URL/webservice..." ) );
        _file_jmenu.addSeparator();
        _file_jmenu.add( _exit_item = new JMenuItem( "Exit" ) );
        customizeJMenuItem( _open_url_item );
        customizeJMenuItem( _exit_item );
        _jmenubar.add( _file_jmenu );
    }

    void buildFontSizeMenu() {
        _font_size_menu = createMenu( FONT_SIZE_MENU_LABEL, getConfiguration() );
        _font_size_menu.add( _super_tiny_fonts_item = new JMenuItem( "Super Tiny Fonts" ) );
        _font_size_menu.add( _tiny_fonts_item = new JMenuItem( "Tiny Fonts" ) );
        _font_size_menu.add( _small_fonts_item = new JMenuItem( "Small Fonts" ) );
        _font_size_menu.add( _medium_fonts_item = new JMenuItem( "Medium Fonts" ) );
        _font_size_menu.add( _large_fonts_item = new JMenuItem( "Large Fonts" ) );
        customizeJMenuItem( _super_tiny_fonts_item );
        customizeJMenuItem( _tiny_fonts_item );
        customizeJMenuItem( _small_fonts_item );
        customizeJMenuItem( _medium_fonts_item );
        customizeJMenuItem( _large_fonts_item );
        _jmenubar.add( _font_size_menu );
    }

    void buildHelpMenu() {
        _help_jmenu = createMenu( "Help", getConfiguration() );
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
        _phyloxml_ref_item.setToolTipText( PHYLOXML_REF_TOOL_TIP );
        _aptx_ref_item.setToolTipText( APTX_REF_TOOL_TIP );
        _jmenubar.add( _help_jmenu );
    }

    void buildToolsMenu() {
        _tools_menu = createMenu( "Tools", getConfiguration() );
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
        _type_menu = createMenu( TYPE_MENU_HEADER, getConfiguration() );
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
        _view_jmenu = createMenu( "View as Text", getConfiguration() );
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
                                  "Please enter the minimum for confidence values to be displayed.\n"
                                          + "[current value: " + getOptions().getMinConfidenceValue() + "]\n",
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

    void close() {
        removeTextFrame();
        if ( _mainpanel != null ) {
            _mainpanel.terminate();
        }
        if ( _contentpane != null ) {
            _contentpane.removeAll();
        }
        setVisible( false );
        dispose();
    }

    void confColor() {
        if ( _mainpanel.getCurrentTreePanel() != null ) {
            _mainpanel.getCurrentTreePanel().confColor();
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
        if ( jmi != null ) {
            jmi.setFont( MainFrame.menu_font );
            if ( !getConfiguration().isUseNativeUI() ) {
                jmi.setBackground( Constants.MENU_BACKGROUND_COLOR_DEFAULT );
                jmi.setForeground( Constants.MENU_TEXT_COLOR_DEFAULT );
            }
            jmi.addActionListener( this );
        }
    }

    void customizeRadioButtonMenuItem( final JRadioButtonMenuItem item, final boolean is_selected ) {
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

    void exceptionOccuredDuringOpenFile( final Exception e ) {
        try {
            _mainpanel.getCurrentTreePanel().setArrowCursor();
        }
        catch ( final Exception ex ) {
            // Do nothing.
        }
        JOptionPane.showMessageDialog( this,
                                       ForesterUtil.wordWrap( e.getLocalizedMessage(), 80 ),
                                       "Error during File|Open",
                                       JOptionPane.ERROR_MESSAGE );
    }

    void exceptionOccuredDuringSaveAs( final Exception e ) {
        try {
            _mainpanel.getCurrentTreePanel().setArrowCursor();
        }
        catch ( final Exception ex ) {
            // Do nothing.
        }
        JOptionPane.showMessageDialog( this, "Exception" + e, "Error during File|SaveAs", JOptionPane.ERROR_MESSAGE );
    }

    boolean GAndSDoHaveMoreThanOneSpeciesInComman( final Phylogeny gene_tree ) {
        if ( ( gene_tree == null ) || gene_tree.isEmpty() ) {
            JOptionPane.showMessageDialog( this,
                                           "Gene tree and species tree have no species in common.",
                                           "Error during SDI",
                                           JOptionPane.ERROR_MESSAGE );
            return false;
        }
        else if ( gene_tree.getNumberOfExternalNodes() < 2 ) {
            JOptionPane.showMessageDialog( this,
                                           "Gene tree and species tree have only one species in common.",
                                           "Error during SDI",
                                           JOptionPane.ERROR_MESSAGE );
            return false;
        }
        else {
            return true;
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

    MainPanel getMainPanel() {
        return _mainpanel;
    }

    Options getOptions() {
        return _options;
    }

    void initializeTypeMenu( final Options options ) {
        setTypeMenuToAllUnselected();
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

    void midpointRoot() {
        if ( _mainpanel.getCurrentTreePanel() != null ) {
            _mainpanel.getCurrentTreePanel().midpointRoot();
        }
    }

    abstract void readPhylogeniesFromURL();

    void readPhylogeniesFromWebservice( final int i ) {
        final UrlTreeReader reader = new UrlTreeReader( this, i );
        new Thread( reader ).start();
    }

    void readPhylogeniesFromWebserviceOldyButGoody( final int i ) {
        URL url = null;
        Phylogeny[] trees = null;
        final WebservicesManager webservices_manager = WebservicesManager.getInstance();
        final PhylogeniesWebserviceClient client = webservices_manager.getAvailablePhylogeniesWebserviceClient( i );
        String identifier = JOptionPane.showInputDialog( this, client.getInstructions() + "\n(Reference: "
                + client.getReference() + ")", client.getDescription(), JOptionPane.QUESTION_MESSAGE );
        if ( ( identifier != null ) && ( identifier.trim().length() > 0 ) ) {
            identifier = identifier.trim();
            if ( client.isQueryInteger() ) {
                int id = -1;
                try {
                    id = Integer.parseInt( identifier );
                }
                catch ( final NumberFormatException e ) {
                    id = -1;
                }
                if ( id < 1 ) {
                    JOptionPane.showMessageDialog( this,
                                                   "Identifier is expected to be a number",
                                                   "Can not open URL",
                                                   JOptionPane.ERROR_MESSAGE );
                    return;
                }
                identifier = id + "";
            }
            try {
                String url_str = client.getUrl();
                url_str = url_str.replaceFirst( PhylogeniesWebserviceClient.QUERY_PLACEHOLDER, identifier );
                url = new URL( url_str );
                PhylogenyParser parser = null;
                switch ( client.getReturnFormat() ) {
                    case TOL_XML_RESPONSE:
                        parser = new TolParser();
                        break;
                    case NEXUS:
                        parser = new NexusPhylogeniesParser();
                        break;
                    case NH:
                        parser = new NHXParser();
                        ( ( NHXParser ) parser ).setTaxonomyExtraction( TAXONOMY_EXTRACTION.NO );
                        ( ( NHXParser ) parser ).setReplaceUnderscores( true );
                        ( ( NHXParser ) parser ).setGuessRootedness( true );
                        break;
                    case NH_EXTRACT_TAXONOMY:
                        parser = new NHXParser();
                        ( ( NHXParser ) parser ).setTaxonomyExtraction( TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
                        ( ( NHXParser ) parser ).setReplaceUnderscores( false );
                        ( ( NHXParser ) parser ).setGuessRootedness( true );
                        break;
                    case PFAM:
                        parser = new NHXParser();
                        ( ( NHXParser ) parser ).setTaxonomyExtraction( TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY );
                        ( ( NHXParser ) parser ).setReplaceUnderscores( false );
                        ( ( NHXParser ) parser ).setGuessRootedness( true );
                        break;
                    case NHX:
                        parser = new NHXParser();
                        ( ( NHXParser ) parser ).setTaxonomyExtraction( TAXONOMY_EXTRACTION.NO );
                        ( ( NHXParser ) parser ).setReplaceUnderscores( false );
                        ( ( NHXParser ) parser ).setGuessRootedness( true );
                        break;
                    case PHYLOXML:
                        parser = new PhyloXmlParser();
                        break;
                    default:
                        throw new IllegalArgumentException( "unknown format: " + client.getReturnFormat() );
                }
                if ( _mainpanel.getCurrentTreePanel() != null ) {
                    _mainpanel.getCurrentTreePanel().setWaitCursor();
                }
                else {
                    _mainpanel.setWaitCursor();
                }
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                trees = factory.create( url.openStream(), parser );
            }
            catch ( final MalformedURLException e ) {
                JOptionPane.showMessageDialog( this,
                                               "Malformed URL: " + url + "\n" + e.getLocalizedMessage(),
                                               "Malformed URL",
                                               JOptionPane.ERROR_MESSAGE );
            }
            catch ( final IOException e ) {
                JOptionPane.showMessageDialog( this,
                                               "Could not read from " + url + "\n" + e.getLocalizedMessage(),
                                               "Failed to read tree from " + client.getName() + " for " + identifier,
                                               JOptionPane.ERROR_MESSAGE );
            }
            catch ( final NumberFormatException e ) {
                JOptionPane.showMessageDialog( this,
                                               "Could not read from " + url + "\n" + e.getLocalizedMessage(),
                                               "Failed to read tree from " + client.getName() + " for " + identifier,
                                               JOptionPane.ERROR_MESSAGE );
            }
            catch ( final Exception e ) {
                e.printStackTrace();
                JOptionPane.showMessageDialog( this,
                                               e.getLocalizedMessage(),
                                               "Unexpected Exception",
                                               JOptionPane.ERROR_MESSAGE );
            }
            finally {
                if ( _mainpanel.getCurrentTreePanel() != null ) {
                    _mainpanel.getCurrentTreePanel().setArrowCursor();
                }
                else {
                    _mainpanel.setArrowCursor();
                }
            }
            if ( ( trees != null ) && ( trees.length > 0 ) ) {
                for( final Phylogeny phylogeny : trees ) {
                    if ( client.getName().equals( WebserviceUtil.TREE_FAM_NAME ) ) {
                        phylogeny.setRerootable( false );
                    }
                    if ( client.getName().equals( WebserviceUtil.PFAM_NAME ) ) {
                        phylogeny.setRerootable( false );
                        ForesterUtil.transferInternalNodeNamesToConfidence( phylogeny );
                    }
                    if ( client.getProcessingInstructions() != null ) {
                        WebserviceUtil.processInstructions( client, phylogeny );
                    }
                    if ( client.getNodeField() != null ) {
                        ForesterUtil.transferNodeNameToField( phylogeny, client.getNodeField() );
                    }
                    phylogeny.setIdentifier( new Identifier( identifier, client.getName() ) );
                    _jmenubar.remove( _help_jmenu );
                    _jmenubar.add( _help_jmenu );
                    _mainpanel.addPhylogenyInNewTab( phylogeny,
                                                     getConfiguration(),
                                                     new File( url.getFile() ).getName(),
                                                     url.toString() );
                    String my_name_for_file = "";
                    if ( !ForesterUtil.isEmpty( phylogeny.getName() ) ) {
                        my_name_for_file = new String( phylogeny.getName() ).replaceAll( " ", "_" );
                    }
                    else if ( phylogeny.getIdentifier() != null ) {
                        final StringBuffer sb = new StringBuffer();
                        if ( !ForesterUtil.isEmpty( phylogeny.getIdentifier().getProvider() ) ) {
                            sb.append( phylogeny.getIdentifier().getProvider() );
                            sb.append( "_" );
                        }
                        sb.append( phylogeny.getIdentifier().getValue() );
                        my_name_for_file = new String( sb.toString().replaceAll( " ", "_" ) );
                    }
                    _mainpanel.getCurrentTreePanel().setTreeFile( new File( my_name_for_file ) );
                    Util.lookAtSomeTreePropertiesForAptxControlSettings( phylogeny,
                                                                         _mainpanel.getControlPanel(),
                                                                         getConfiguration() );
                    _mainpanel.getControlPanel().showWhole();
                }
            }
            _contentpane.repaint();
        }
        activateSaveAllIfNeeded();
        System.gc();
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

    void setOptions( final Options options ) {
        _options = options;
    }

    void setSelectedTypeInTypeMenu( final PHYLOGENY_GRAPHICS_TYPE type ) {
        setTypeMenuToAllUnselected();
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

    void setTypeMenuToAllUnselected() {
        _convex_type_cbmi.setSelected( false );
        _curved_type_cbmi.setSelected( false );
        _euro_type_cbmi.setSelected( false );
        _rounded_type_cbmi.setSelected( false );
        _triangular_type_cbmi.setSelected( false );
        _rectangular_type_cbmi.setSelected( false );
        _unrooted_type_cbmi.setSelected( false );
        _circular_type_cbmi.setSelected( false );
    }

    void showWhole() {
        _mainpanel.getControlPanel().showWhole();
    }

    void switchColors() {
        final TreeColorSet colorset = _mainpanel.getTreeColorSet();
        final ColorSchemeChooser csc = new ColorSchemeChooser( getMainPanel(), colorset );
        csc.setVisible( true );
    }

    void taxColor() {
        if ( _mainpanel.getCurrentTreePanel() != null ) {
            _mainpanel.getCurrentTreePanel().taxColor();
        }
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
            updateScreenTextAntialias( getMainPanel().getTreePanels() );
            if ( getCurrentTreePanel().getControlPanel().getDynamicallyHideData() != null ) {
                if ( new_type == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
                    getCurrentTreePanel().getControlPanel().getDynamicallyHideData().setEnabled( false );
                }
                else {
                    getCurrentTreePanel().getControlPanel().getDynamicallyHideData().setEnabled( true );
                }
            }
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
        options.setPrintUsingActualSize( ( _print_using_actual_size_cbmi != null )
                && ( _print_using_actual_size_cbmi.isSelected() ) );
        options.setGraphicsExportUsingActualSize( ( _graphics_export_using_actual_size_cbmi != null )
                && ( _graphics_export_using_actual_size_cbmi.isSelected() ) );
        options.setAntialiasPrint( ( _antialias_print_cbmi != null ) && _antialias_print_cbmi.isSelected() );
        options.setPrintBlackAndWhite( ( _print_black_and_white_cbmi != null )
                && _print_black_and_white_cbmi.isSelected() );
        options
                .setInternalNumberAreConfidenceForNhParsing( ( _internal_number_are_confidence_for_nh_parsing_cbmi != null )
                        && _internal_number_are_confidence_for_nh_parsing_cbmi.isSelected() );
        options.setNhParsingIgnoreQuotes( ( _ignore_quotes_in_nh_parsing_cbmi != null )
                && _ignore_quotes_in_nh_parsing_cbmi.isSelected() );
        options.setExtractPfamTaxonomyCodesInNhParsing( ( _extract_pfam_style_tax_codes_cbmi != null )
                && _extract_pfam_style_tax_codes_cbmi.isSelected() );
        options.setReplaceUnderscoresInNhParsing( ( _replace_underscores_cbmi != null )
                && _replace_underscores_cbmi.isSelected() );
        options.setMatchWholeTermsOnly( ( _search_whole_words_only_cbmi != null )
                && _search_whole_words_only_cbmi.isSelected() );
        options.setInverseSearchResult( ( _inverse_search_result_cbmi != null )
                && _inverse_search_result_cbmi.isSelected() );
        if ( _graphics_export_visible_only_cbmi != null ) {
            options.setGraphicsExportVisibleOnly( _graphics_export_visible_only_cbmi.isSelected() );
            if ( _graphics_export_visible_only_cbmi.isSelected() && ( _graphics_export_using_actual_size_cbmi != null ) ) {
                _graphics_export_using_actual_size_cbmi.setSelected( true );
                _graphics_export_using_actual_size_cbmi.setEnabled( false );
            }
            else {
                _graphics_export_using_actual_size_cbmi.setEnabled( true );
            }
        }
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
        if ( ( _mainpanel.getCurrentPhylogeny() == null ) || _mainpanel.getCurrentPhylogeny().isEmpty()
                || ( _mainpanel.getCurrentPhylogeny().getNumberOfExternalNodes() > 10000 ) ) {
            return;
        }
        _textframe = TextFrame.instantiate( _mainpanel.getCurrentPhylogeny().toNexus() );
    }

    void viewAsNH() {
        removeTextFrame();
        if ( ( _mainpanel.getCurrentPhylogeny() == null ) || _mainpanel.getCurrentPhylogeny().isEmpty()
                || ( _mainpanel.getCurrentPhylogeny().getNumberOfExternalNodes() > 10000 ) ) {
            return;
        }
        _textframe = TextFrame.instantiate( _mainpanel.getCurrentPhylogeny().toNewHampshire( false ) );
    }

    void viewAsNHX() {
        removeTextFrame();
        if ( ( _mainpanel.getCurrentPhylogeny() == null ) || _mainpanel.getCurrentPhylogeny().isEmpty()
                || ( _mainpanel.getCurrentPhylogeny().getNumberOfExternalNodes() > 10000 ) ) {
            return;
        }
        _textframe = TextFrame.instantiate( _mainpanel.getCurrentPhylogeny().toNewHampshireX() );
    }

    void viewAsXML() {
        removeTextFrame();
        if ( ( _mainpanel.getCurrentPhylogeny() != null ) && !_mainpanel.getCurrentPhylogeny().isEmpty()
                && ( _mainpanel.getCurrentPhylogeny().getNumberOfExternalNodes() <= 10000 ) ) {
            _textframe = TextFrame.instantiate( _mainpanel.getCurrentPhylogeny().toPhyloXML( 0 ) );
        }
    }

    /**
     * Display the about box.
     */
    static void about() {
        final StringBuffer about = new StringBuffer( "Archaeopteryx\nVersion " + Constants.VERSION + "\n" );
        about.append( "Copyright (C) 2007-2009 Christian M Zmasek\n" );
        about.append( "Copyright (C) 2007-2009 Ethalinda KS Cannon\n" );
        about.append( "All Rights Reserved\n" );
        about.append( "License: GNU Lesser General Public License (LGPL)\n" );
        about.append( "Last modified: " + Constants.PRG_DATE + "\n" );
        about.append( "phyloXML version : " + ForesterConstants.PHYLO_XML_VERSION + "\n" );
        about.append( "phyloXML location: " + ForesterConstants.PHYLO_XML_LOCATION + "\n" );
        if ( !ForesterUtil.isEmpty( ForesterUtil.JAVA_VERSION ) && !ForesterUtil.isEmpty( ForesterUtil.JAVA_VENDOR ) ) {
            about.append( "[your Java version: " + ForesterUtil.JAVA_VERSION + " " + ForesterUtil.JAVA_VENDOR + "]\n" );
        }
        if ( !ForesterUtil.isEmpty( ForesterUtil.OS_NAME ) && !ForesterUtil.isEmpty( ForesterUtil.OS_ARCH )
                && !ForesterUtil.isEmpty( ForesterUtil.OS_VERSION ) ) {
            about.append( "[your OS: " + ForesterUtil.OS_NAME + " " + ForesterUtil.OS_ARCH + " "
                    + ForesterUtil.OS_VERSION + "]\n" );
        }
        final Runtime rt = java.lang.Runtime.getRuntime();
        final long free_memory = rt.freeMemory() / 1000000;
        final long total_memory = rt.totalMemory() / 1000000;
        about.append( "[free memory: " + free_memory + "MB, total memory: " + total_memory + "MB]\n" );
        about.append( "[locale: " + Locale.getDefault() + "]\n" );
        about.append( "References:\n" );
        about.append( Constants.PHYLOXML_REFERENCE_SHORT + "\n" );
        about.append( "Zmasek CM and Eddy SR (2001), Bioinformatics, 17, 383\n" );
        about.append( "For more information & download:\n" );
        about.append( Constants.APTX_WEB_SITE + "\n" );
        about.append( "Comments: " + Constants.AUTHOR_EMAIL );
        JOptionPane.showMessageDialog( null, about, Constants.PRG_NAME, JOptionPane.PLAIN_MESSAGE );
    }

    static String createCurrentFontDesc( final TreeFontSet tree_font_set ) {
        return tree_font_set.getLargeFont().getFamily() + " " + tree_font_set.getLargeFont().getSize();
    }

    static JMenu createMenu( final String title, final Configuration conf ) {
        final JMenu jmenu = new JMenu( title );
        if ( !conf.isUseNativeUI() ) {
            jmenu.setFont( MainFrame.menu_font );
            jmenu.setBackground( Constants.MENU_BACKGROUND_COLOR_DEFAULT );
            jmenu.setForeground( Constants.MENU_TEXT_COLOR_DEFAULT );
        }
        return jmenu;
    }

    static JMenuItem customizeMenuItemAsLabel( final JMenuItem label, final Configuration configuration ) {
        label.setFont( MainFrame.menu_font.deriveFont( Font.BOLD ) );
        if ( !configuration.isUseNativeUI() ) {
            label.setBackground( Constants.MENU_LABEL_BACKGROUND_COLOR_DEFAULT );
            label.setForeground( Constants.MENU_LABEL_TEXT_COLOR_DEFAULT );
            label.setOpaque( true );
        }
        label.setSelected( false );
        label.setEnabled( false );
        return label;
    }

    static void cycleOverview( final Options op, final TreePanel tree_panel ) {
        switch ( op.getOvPlacement() ) {
            case LOWER_LEFT:
                op.setOvPlacement( Options.OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT );
                break;
            case LOWER_RIGHT:
                op.setOvPlacement( Options.OVERVIEW_PLACEMENT_TYPE.LOWER_LEFT );
                break;
            case UPPER_LEFT:
                op.setOvPlacement( Options.OVERVIEW_PLACEMENT_TYPE.UPPER_RIGHT );
                break;
            case UPPER_RIGHT:
                op.setOvPlacement( Options.OVERVIEW_PLACEMENT_TYPE.LOWER_RIGHT );
                break;
            default:
                throw new IllegalStateException( "unknown placement: " + op.getOvPlacement() );
        }
        tree_panel.updateOvSettings();
    }

    static void help( final Map<String, WebLink> weblinks ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( "Display options\n" );
        sb.append( "-------------------\n" );
        sb.append( "Use the checkboxes to select types of information to display on the tree.\n\n" );
        sb.append( "Clickable tree nodes\n" );
        sb.append( "--------------------\n" );
        sb.append( "Tree nodes can be clicked, the action is determined by the 'click on node to' menu\n" );
        sb.append( "or by right clicking:\n" );
        sb.append( "o  Display Node Data -- display information for a node\n" );
        sb.append( "o  Collapse/Uncollapse -- collapse and uncollapse subtree from clicked node\n" );
        sb.append( "o  Root/Reroot -- change tree root to clicked node\n" );
        sb.append( "o  Sub/Super Tree -- toggle between subtree from clicked node and whole tree\n" );
        sb.append( "o  Swap Descendants -- switch descendant on either side of clicked node\n" );
        sb.append( "o  Colorize Subtree -- color a subtree\n" );
        sb.append( "o  Open Sequence Web -- launch a web browser to display sequence information\n" );
        sb.append( "o  Open Taxonomy Web -- launch a web browser to display taxonomy information\n" );
        sb.append( "-  there may be additional choices depending on this particular setup\n\n" );
        sb.append( "Right clicking on a node always displays the information of a node.\n\n" );
        sb.append( "Zooming\n" );
        sb.append( "---------\n" );
        sb.append( "The mouse wheel and the plus and minus keys control zooming.\n" );
        sb.append( "Mouse wheel+Ctrl changes the text size.\n" );
        sb.append( "Mouse wheel+Shift controls zooming in vertical direction only.\n" );
        sb.append( "Use the buttons on the control panel to zoom the tree in and out, horizontally or vertically.\n" );
        sb
                .append( "The entire tree can be fitted into the window by clicking the \"F\" button, or by pressing F, Delete, or Home.\n" );
        sb.append( "The up, down, left, and right keys can be used to move the visible part (if zoomed in).\n" );
        sb.append( "Up, down, left, and right+Shift can be used to control zooming horizontally and vertically.\n" );
        sb.append( "Plus and minus keys+Ctrl change the text size; F+Ctrl, Delete+Ctrl, or Home+Ctrl resets it.\n\n" );
        sb.append( "Quick tree manipulation:\n" );
        sb.append( "------------------------\n" );
        sb.append( "Order Subtrees -- order the tree by branch length\n" );
        sb.append( "Uncollapse All -- uncollapse any and all collapsed branches\n\n" );
        sb.append( "Memory problems (Java heap space error)\n" );
        sb.append( "---------------------------------------\n" );
        sb.append( "Since the Java default memory allocation is quite small, it might by necessary (for trees\n" );
        sb
                .append( "with more than approximately 5000 external nodes) to increase the memory which Java can use, with\n" );
        sb.append( "the '-Xmx' Java command line option. For example:\n" );
        sb.append( "java -Xms32m -Xmx256m -cp path\\to\\forester.jar org.forester.archaeopteryx.Archaeopteryx\n\n" );
        if ( ( weblinks != null ) && ( weblinks.size() > 0 ) ) {
            sb.append( "Active web links\n" );
            sb.append( "--------------------\n" );
            for( final String key : weblinks.keySet() ) {
                sb.append( " " + weblinks.get( key ).toString() + "\n" );
            }
        }
        // + "General remarks\n"
        // + "---------------\n"
        // +
        // "o  The application version permits copying to the clipboard \n"
        // +
        // "    in the \"View\"|\"View as ...\" frame (either by control-c or button press).\n"
        // +
        // "o  Changes made to a subtree affect this subtree and its subtrees,\n"
        // + "    but not any of its parent tree(s).\n"
        // +
        // "o  Archaeopteryx tries to detect whether the numerical values in a NH tree\n"
        // +
        // "    are likely to be bootstrap values instead of branch length values.\n\n"
        // +
        // " Remarks regarding SDI (Speciation Duplication Inference):\n"
        // +
        // "o  Each external node of the gene tree (in display) needs to be associated with\n"
        // +
        // "    a species: either directly through the \"Species\" field, or the species\n"
        // +
        // "    is part of the sequence name in the form \"XXXX_SPECIES\"\n"
        // +
        // "    (e.g. \"ACON_DROME\" or \"ACON_DROME/123-4489\" which is also acceptable).\n"
        // +
        // "o  A species tree for each species of the gene tree needs to be loaded with\n"
        // +
        // "   \"SDI\"|\"Load species tree\" prior the SDI execution.\n"
        // +
        // "o  !External nodes of the gene tree associated with species not present in\n"
        // +
        // "    the species tree are REMOVED prior to SDI execution!\n"
        // +
        // "o  Both the gene tree and the species tree must be completely binary.\n"
        // +
        // "o  Duplications and speciations are a function of the position of the root.\n"
        // +
        // "    Hence, after each manual \"Root/Reroot\"ing some duplications will be\n"
        // + "    incorrect and need to be inferred again\n"
        // +
        // "    with: \"SDI\"|\"SDI (Speciation Duplication Inference)\n\n"
        sb.append( "\n" );
        sb.append( "phyloXML\n" );
        sb.append( "-------------------\n" );
        sb.append( "Reference: " + Constants.PHYLOXML_REFERENCE + "\n" );
        sb.append( "Website: " + Constants.PHYLOXML_WEB_SITE + "\n" );
        sb.append( "Version: " + ForesterConstants.PHYLO_XML_VERSION + "\n" );
        sb.append( "\n" );
        sb.append( "For more information: http://www.phylosoft.org/archaeopteryx/\n" );
        sb.append( "Email: " + Constants.AUTHOR_EMAIL + "\n\n" );
        TextFrame.instantiate( sb.toString() );
    }

    static void setOvPlacementColorChooseMenuItem( final JMenuItem mi, final TreePanel tree_panel ) {
        if ( ( tree_panel != null ) && ( tree_panel.getTreeColorSet() != null ) ) {
            mi.setText( "Overview Placement... (current: " + tree_panel.getOptions().getOvPlacement() + ")" );
        }
        else {
            mi.setText( "Overview Placement..." );
        }
    }

    static void setTextColorChooseMenuItem( final JMenuItem mi, final TreePanel tree_panel ) {
        if ( ( tree_panel != null ) && ( tree_panel.getTreeColorSet() != null ) ) {
            mi.setText( "Select Colors... (current: " + tree_panel.getTreeColorSet().getCurrentColorSchemeName() + ")" );
        }
        else {
            mi.setText( "Select Colors..." );
        }
    }

    static void setTextForFontChooserMenuItem( final JMenuItem mi, final String font_desc ) {
        mi.setText( "Select Font... (current: " + font_desc + ")" );
    }

    static void setTextMinSupportMenuItem( final JMenuItem mi, final Options options, final TreePanel current_tree_panel ) {
        if ( ( current_tree_panel == null ) || ( current_tree_panel.getPhylogeny() == null ) ) {
            mi.setEnabled( true );
        }
        else if ( ForesterUtil.isHasAtLeastOneBranchWithSupportValues( current_tree_panel.getPhylogeny() ) ) {
            mi.setEnabled( true );
        }
        else {
            mi.setEnabled( false );
        }
        mi.setText( "Enter Min Confidence Value... (current: " + options.getMinConfidenceValue() + ")" );
    }

    static void updateOptionsMenuDependingOnPhylogenyType( final MainPanel main_panel,
                                                           final JCheckBoxMenuItem scale,
                                                           final JCheckBoxMenuItem branch_lengths,
                                                           final JRadioButtonMenuItem non_lined_up,
                                                           final JRadioButtonMenuItem uniform_clado,
                                                           final JRadioButtonMenuItem nonuniform_clado,
                                                           final JCheckBoxMenuItem label_direction_cbmi ) {
        final TreePanel tree_panel = main_panel.getCurrentTreePanel();
        final ControlPanel control = main_panel.getControlPanel();
        final Options options = main_panel.getOptions();
        scale.setSelected( options.isShowScale() );
        branch_lengths.setSelected( options.isShowBranchLengthValues() );
        // non_lined_up.setSelected( options.isNonLinedUpCladogram() );
        if ( ( tree_panel != null ) && ( !tree_panel.isPhyHasBranchLengths() ) ) {
            scale.setSelected( false );
            scale.setEnabled( false );
            branch_lengths.setSelected( false );
            branch_lengths.setEnabled( false );
        }
        else if ( ( tree_panel != null ) && !control.isDrawPhylogram() ) {
            scale.setSelected( false );
            scale.setEnabled( false );
            branch_lengths.setEnabled( true );
        }
        else {
            scale.setEnabled( true );
            branch_lengths.setEnabled( true );
        }
        if ( ( tree_panel != null )
                && ( ( tree_panel.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.ROUNDED )
                        && ( tree_panel.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) && ( tree_panel
                        .getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) ) ) {
            branch_lengths.setSelected( false );
            branch_lengths.setEnabled( false );
        }
        if ( tree_panel != null ) {
            if ( ( tree_panel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR )
                    || ( tree_panel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) ) {
                non_lined_up.setEnabled( false );
                uniform_clado.setEnabled( false );
                nonuniform_clado.setEnabled( false );
            }
            else {
                non_lined_up.setEnabled( true );
                uniform_clado.setEnabled( true );
                nonuniform_clado.setEnabled( true );
            }
        }
        else {
            if ( ( tree_panel != null )
                    && ( ( tree_panel.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) && ( tree_panel
                            .getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) ) ) {
                branch_lengths.setSelected( false );
                branch_lengths.setEnabled( false );
            }
            if ( ( tree_panel != null )
                    && ( ( tree_panel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) || ( tree_panel
                            .getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) ) ) {
                non_lined_up.setEnabled( false );
            }
            else {
                // non_lined_up.setSelected( options.isNonLinedUpCladogram() );
                non_lined_up.setEnabled( true );
            }
        }
        label_direction_cbmi.setEnabled( true );
        if ( tree_panel != null ) {
            if ( ( tree_panel.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                    && ( tree_panel.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                label_direction_cbmi.setEnabled( false );
            }
            if ( tree_panel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                scale.setSelected( false );
                scale.setEnabled( false );
            }
        }
    }

    static void updateScreenTextAntialias( final List<TreePanel> treepanels ) {
        for( final TreePanel tree_panel : treepanels ) {
            tree_panel.setTextAntialias();
        }
    }
}
