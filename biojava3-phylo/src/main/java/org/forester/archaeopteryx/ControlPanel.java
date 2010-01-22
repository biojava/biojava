// $Id: ControlPanel.java,v 1.49 2009/11/18 19:07:28 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.JTextField;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

final class ControlPanel extends JPanel implements ActionListener {

    private static final String  RETURN_TO_SUPER_TREE_TEXT = "Back to Super Tree";
    final static Font            jcb_font                  = new Font( Configuration.getDefaultFontFamilyName(),
                                                                       Font.PLAIN,
                                                                       9 );
    final static Font            js_font                   = new Font( Configuration.getDefaultFontFamilyName(),
                                                                       Font.PLAIN,
                                                                       9 );
    final static Font            jcb_bold_font             = new Font( Configuration.getDefaultFontFamilyName(),
                                                                       Font.BOLD,
                                                                       9 );
    private static Color         background_color;
    private static Color         jcb_text_color;
    private static Color         jcb_background_color;
    private static Color         button_text_color;
    private static Color         button_background_color;
    private static Color         button_border_color;
    private static final long    serialVersionUID          = -8463483932821545633L;
    private final MainPanel      _mainpanel;
    // The settings from the conf file
    private final Configuration  _configuration;
    // Tree checkboxes
    private JCheckBox            _display_internal_data;
    private JCheckBox            _show_node_names;
    private JCheckBox            _show_taxo_code;
    private JCheckBox            _write_confidence;
    private JCheckBox            _show_events;
    private JCheckBox            _color_acc_species;
    private JCheckBox            _color_branches_cb;
    private JCheckBox            _width_branches;
    private JCheckBox            _show_domain_architectures;
    private JCheckBox            _show_annotation;
    private JCheckBox            _show_binary_characters;
    private JCheckBox            _show_binary_character_counts;
    private JCheckBox            _show_gene_names;
    private JCheckBox            _show_gene_symbols;
    private JCheckBox            _show_sequence_acc;
    private JCheckBox            _node_desc_popup_cb;
    private JCheckBox            _dynamically_hide_data;
    private JCheckBox            _show_taxo_names;
    private JCheckBox            _color_according_to_annotation;
    private JCheckBox            _display_as_phylogram_cb;
    private JLabel               _click_to_label;
    private JLabel               _zoom_label;
    private JLabel               _domain_display_label;
    private JComboBox            _click_to_combobox;
    private Map<Integer, String> _all_click_to_names;
    private List<String>         _click_to_names;
    // Indices for the click-to options in the combo box
    private int                  _show_data_item;
    private int                  _collapse_cb_item;
    private int                  _reroot_cb_item;
    private int                  _swap_cb_item;
    private int                  _subtree_cb_item;
    private int                  _color_subtree_cb_item;
    private int                  _open_seq_web_item;
    private int                  _open_tax_web_item;
    private int                  _cut_subtree_item;
    private int                  _copy_subtree_item;
    private int                  _delete_node_or_subtree_item;
    private int                  _paste_subtree_item;
    private int                  _add_new_node_item;
    private int                  _edit_node_data_item;
    private int                  _blast_item;
    // zooming and quick tree manipulation buttons:
    private JButton              _zoom_in_x;
    private JButton              _zoom_in_y;
    private JButton              _zoom_out_x;
    private JButton              _zoom_out_y;
    private JButton              _show_whole;
    private JButton              _order;
    private JButton              _uncollapse_all;
    private JButton              _zoom_in_domain_structure;
    private JButton              _zoom_out_domain_structure;
    private JButton              _decr_domain_structure_evalue_thr;
    private JButton              _incr_domain_structure_evalue_thr;
    private JButton              _return_to_super_tree;
    private JTextField           _domain_structure_evalue_thr_tf;
    private JTextField           _search_tf;
    private boolean              _order_of_appearance;
    private boolean              _color_branches;
    private NodeClickAction      _action_when_node_clicked;
    private List<Boolean>        _draw_phylogram;
    private Map<String, Color>   _annotation_colors;
    private Map<String, Color>   _species_colors;
    private JButton              _search_reset_button;
    private JLabel               _search_found_label;

    ControlPanel( final MainPanel ap, final Configuration config_settings ) {
        init();
        _mainpanel = ap;
        _configuration = config_settings;
        setDefaultColors();
        if ( !_configuration.isUseNativeUI() ) {
            setBackground( ControlPanel.background_color );
            setBorder( BorderFactory.createRaisedBevelBorder() );
        }
        setLayout( new GridLayout( 0, 1, 2, 2 ) );
        _order_of_appearance = true;
        setupControls();
        //Dimension dim =new Dimension();
        //dim.height =200;
        //dim.width =200;
        //setPreferredSize( dim);
    }

    /**
     * Handle an action.
     */
    public void actionPerformed( final ActionEvent e ) {
        try {
            final TreePanel tp = _mainpanel.getCurrentTreePanel();
            if ( tp == null ) {
                return;
            }
            if ( e.getSource() == _click_to_combobox ) {
                setClickToAction( _click_to_combobox.getSelectedIndex() );
                getCurrentTreePanel().repaint();
            }
            else if ( e.getSource() == _show_binary_characters ) {
                if ( ( _show_binary_character_counts != null ) && _show_binary_characters.isSelected() ) {
                    _show_binary_character_counts.setSelected( false );
                }
                displayedPhylogenyMightHaveChanged( true );
            }
            else if ( e.getSource() == _show_binary_character_counts ) {
                if ( ( _show_binary_characters != null ) && _show_binary_character_counts.isSelected() ) {
                    _show_binary_characters.setSelected( false );
                }
                displayedPhylogenyMightHaveChanged( true );
            }
            else if ( e.getSource() == _color_according_to_annotation ) {
                if ( ( _show_annotation != null ) && _color_according_to_annotation.isSelected() ) {
                    _show_annotation.setSelected( true );
                }
                displayedPhylogenyMightHaveChanged( false );
            }
            else if ( e.getSource() == _show_annotation ) {
                if ( ( _color_according_to_annotation != null ) && !_show_annotation.isSelected() ) {
                    _color_according_to_annotation.setSelected( false );
                }
                displayedPhylogenyMightHaveChanged( false );
            }
            else if ( ( tp != null ) && ( tp.getPhylogeny() != null ) ) {
                if ( e.getSource() == getDisplayAsPhylogramCb() ) {
                    setDrawPhylogram( getDisplayAsPhylogramCb().isSelected() );
                    showWhole();
                }
                // Zoom buttons
                else if ( e.getSource() == _zoom_in_x ) {
                    zoomInX( Constants.BUTTON_ZOOM_IN_FACTOR, Constants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR );
                    displayedPhylogenyMightHaveChanged( false );
                }
                else if ( e.getSource() == _zoom_in_y ) {
                    zoomInY( Constants.BUTTON_ZOOM_IN_FACTOR );
                    displayedPhylogenyMightHaveChanged( false );
                }
                else if ( e.getSource() == _zoom_out_x ) {
                    zoomOutX( Constants.BUTTON_ZOOM_OUT_FACTOR, Constants.BUTTON_ZOOM_OUT_X_CORRECTION_FACTOR );
                    displayedPhylogenyMightHaveChanged( false );
                }
                else if ( e.getSource() == _zoom_out_y ) {
                    zoomOutY( Constants.BUTTON_ZOOM_OUT_FACTOR );
                    displayedPhylogenyMightHaveChanged( false );
                }
                else if ( e.getSource() == _show_whole ) {
                    showWhole();
                }
                else if ( e.getSource() == _return_to_super_tree ) {
                    _mainpanel.getCurrentTreePanel().superTree();
                    showWhole();
                }
                else if ( e.getSource() == _order ) {
                    tp.getPhylogeny().orderAppearance( _order_of_appearance );
                    _order_of_appearance = !_order_of_appearance;
                    displayedPhylogenyMightHaveChanged( false );
                }
                else if ( e.getSource() == _uncollapse_all ) {
                    uncollapseAll( tp );
                    displayedPhylogenyMightHaveChanged( false );
                }
                else if ( e.getSource() == _zoom_in_domain_structure ) {
                    _mainpanel.getCurrentTreePanel().zoomInDomainStructure();
                    displayedPhylogenyMightHaveChanged( true );
                }
                else if ( e.getSource() == _zoom_out_domain_structure ) {
                    _mainpanel.getCurrentTreePanel().zoomOutDomainStructure();
                    displayedPhylogenyMightHaveChanged( true );
                }
                else if ( e.getSource() == _decr_domain_structure_evalue_thr ) {
                    _mainpanel.getCurrentTreePanel().decreaseDomainStructureEvalueThreshold();
                    displayedPhylogenyMightHaveChanged( true );
                }
                else if ( e.getSource() == _incr_domain_structure_evalue_thr ) {
                    _mainpanel.getCurrentTreePanel().increaseDomainStructureEvalueThreshold();
                    displayedPhylogenyMightHaveChanged( true );
                }
                else if ( e.getSource() == _search_tf ) {
                    search();
                    displayedPhylogenyMightHaveChanged( true );
                }
                else {
                    displayedPhylogenyMightHaveChanged( true );
                }
            }
            tp.requestFocus();
            tp.requestFocusInWindow();
            tp.requestFocus();
        }
        catch ( final Exception ex ) {
            Util.unexpectedException( ex );
        }
        catch ( final Error err ) {
            Util.unexpectedError( err );
        }
    }

    void activateButtonToReturnToSuperTree( int index ) {
        --index;
        if ( index > 0 ) {
            _return_to_super_tree.setText( RETURN_TO_SUPER_TREE_TEXT + " " + index );
        }
        else {
            _return_to_super_tree.setText( RETURN_TO_SUPER_TREE_TEXT );
        }
        _return_to_super_tree.setForeground( Constants.BUTTON_TEXT_COLOR_ON_DEFAULT );
        _return_to_super_tree.setEnabled( true );
    }

    /**
     * Add zoom and quick edit buttons. (Last modified 8/9/04)
     */
    void addButtons() {
        final JLabel spacer = new JLabel( "" );
        spacer.setOpaque( false );
        add( spacer );
        final JPanel x_panel = new JPanel( new GridLayout( 1, 1, 0, 0 ) );
        final JPanel y_panel = new JPanel( new GridLayout( 1, 3, 0, 0 ) );
        final JPanel z_panel = new JPanel( new GridLayout( 1, 1, 0, 0 ) );
        if ( !getConfiguration().isUseNativeUI() ) {
            x_panel.setBackground( getBackground() );
            y_panel.setBackground( getBackground() );
            z_panel.setBackground( getBackground() );
        }
        add( _zoom_label = new JLabel( "Zoom:" ) );
        customizeLabel( _zoom_label, getConfiguration() );
        add( x_panel );
        add( y_panel );
        add( z_panel );
        if ( getConfiguration().isUseNativeUI() ) {
            _zoom_in_x = new JButton( "+" );
            _zoom_out_x = new JButton( "-" );
        }
        else {
            _zoom_in_x = new JButton( "X+" );
            _zoom_out_x = new JButton( "X-" );
        }
        _zoom_in_y = new JButton( "Y+" );
        _zoom_out_y = new JButton( "Y-" );
        _show_whole = new JButton( "F" );
        _show_whole.setToolTipText( "To fit the complete phylogeny to the current display size [Backspace]" );
        _zoom_in_x.setToolTipText( "To zoom in horizontally [Shift+Right]" );
        _zoom_in_y.setToolTipText( "To zoom in vertically [Shift+Up]" );
        _zoom_out_x.setToolTipText( "To zoom out horizontally [Shift+Left]" );
        _zoom_out_y.setToolTipText( "To zoom out vertically [Shift+Down]" );
        if ( getConfiguration().isUseNativeUI() && Util.isMac() ) {
            _zoom_out_x.setPreferredSize( new Dimension( 55, 10 ) );
            _zoom_in_x.setPreferredSize( new Dimension( 55, 10 ) );
        }
        else {
            _zoom_out_x.setPreferredSize( new Dimension( 10, 10 ) );
            _zoom_in_x.setPreferredSize( new Dimension( 10, 10 ) );
        }
        _zoom_out_y.setPreferredSize( new Dimension( 10, 10 ) );
        _zoom_in_y.setPreferredSize( new Dimension( 10, 10 ) );
        _show_whole.setPreferredSize( new Dimension( 10, 10 ) );
        _return_to_super_tree = new JButton( RETURN_TO_SUPER_TREE_TEXT );
        _return_to_super_tree.setEnabled( false );
        _order = new JButton( "Order Subtrees" );
        _uncollapse_all = new JButton( "Uncollapse All" );
        addJButton( _zoom_in_y, x_panel );
        addJButton( _zoom_out_x, y_panel );
        addJButton( _show_whole, y_panel );
        addJButton( _zoom_in_x, y_panel );
        addJButton( _zoom_out_y, z_panel );
        if ( getConfiguration().doDisplayOption( Configuration.show_domain_architectures ) ) {
            setUpControlsForDomainStrucures();
        }
        final JLabel spacer2 = new JLabel( "" );
        add( spacer2 );
        addJButton( _return_to_super_tree, this );
        addJButton( _order, this );
        addJButton( _uncollapse_all, this );
        final JLabel spacer3 = new JLabel( "" );
        add( spacer3 );
        setVisibilityOfDomainStrucureControls();
    }

    void addCheckbox( final int which, final String title ) {
        final JPanel ch_panel = new JPanel( new BorderLayout( 0, 0 ) );
        switch ( which ) {
            case Configuration.display_as_phylogram:
                _display_as_phylogram_cb = new JCheckBox( title );
                getDisplayAsPhylogramCb().setToolTipText( "To switch between phylogram and cladogram display" );
                addJCheckBox( getDisplayAsPhylogramCb(), ch_panel );
                add( ch_panel );
                break;
            case Configuration.display_internal_data:
                _display_internal_data = new JCheckBox( title );
                _display_internal_data.setToolTipText( "To allow or disallow display of internal labels" );
                addJCheckBox( _display_internal_data, ch_panel );
                add( ch_panel );
                break;
            case Configuration.color_according_to_species:
                _color_acc_species = new JCheckBox( title );
                _color_acc_species
                        .setToolTipText( "To colorize taxonomy and sequence labels as a function of taxonomy" );
                addJCheckBox( _color_acc_species, ch_panel );
                add( ch_panel );
                break;
            case Configuration.color_according_to_annotation:
                _color_according_to_annotation = new JCheckBox( title );
                _color_according_to_annotation
                        .setToolTipText( "To colorize sequence annotation labels as a function of sequence annotation" );
                addJCheckBox( _color_according_to_annotation, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_node_names:
                _show_node_names = new JCheckBox( title );
                addJCheckBox( _show_node_names, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_tax_code:
                _show_taxo_code = new JCheckBox( title );
                addJCheckBox( _show_taxo_code, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_taxonomy_names:
                _show_taxo_names = new JCheckBox( title );
                addJCheckBox( _show_taxo_names, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_binary_characters:
                _show_binary_characters = new JCheckBox( title );
                addJCheckBox( _show_binary_characters, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_annotation:
                _show_annotation = new JCheckBox( title );
                addJCheckBox( _show_annotation, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_binary_character_counts:
                _show_binary_character_counts = new JCheckBox( title );
                addJCheckBox( _show_binary_character_counts, ch_panel );
                add( ch_panel );
                break;
            case Configuration.write_confidence_values:
                _write_confidence = new JCheckBox( title );
                addJCheckBox( getWriteConfidenceCb(), ch_panel );
                add( ch_panel );
                break;
            case Configuration.write_events:
                _show_events = new JCheckBox( title );
                addJCheckBox( getShowEventsCb(), ch_panel );
                add( ch_panel );
                break;
            case Configuration.color_branches:
                _color_branches_cb = new JCheckBox( title );
                getColorBranchesCb().setToolTipText( "To use branch color values, if present" );
                addJCheckBox( getColorBranchesCb(), ch_panel );
                add( ch_panel );
                break;
            case Configuration.width_branches:
                _width_branches = new JCheckBox( title );
                _width_branches.setToolTipText( "To use branch width values, if present" );
                addJCheckBox( _width_branches, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_domain_architectures:
                _show_domain_architectures = new JCheckBox( title );
                addJCheckBox( _show_domain_architectures, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_gene_names:
                _show_gene_names = new JCheckBox( title );
                addJCheckBox( _show_gene_names, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_gene_symbols:
                _show_gene_symbols = new JCheckBox( title );
                addJCheckBox( _show_gene_symbols, ch_panel );
                add( ch_panel );
                break;
            case Configuration.show_sequence_acc:
                _show_sequence_acc = new JCheckBox( title );
                addJCheckBox( _show_sequence_acc, ch_panel );
                add( ch_panel );
                break;
            case Configuration.dynamically_hide_data:
                _dynamically_hide_data = new JCheckBox( title );
                getDynamicallyHideData().setToolTipText( "To hide labels depending on likely visibility" );
                addJCheckBox( getDynamicallyHideData(), ch_panel );
                add( ch_panel );
                break;
            case Configuration.node_data_popup:
                _node_desc_popup_cb = new JCheckBox( title );
                getNodeDescPopupCb().setToolTipText( "To enable mouse rollover display of basic node data" );
                addJCheckBox( getNodeDescPopupCb(), ch_panel );
                add( ch_panel );
                break;
            default:
                throw new IllegalStateException( "unknown checkbox: " + which );
        }
    }// addCheckbox

    private void addClickToOption( final int which, final String title ) {
        _click_to_combobox.addItem( title );
        _click_to_names.add( title );
        _all_click_to_names.put( new Integer( which ), title );
        if ( !_configuration.isUseNativeUI() ) {
            _click_to_combobox.setBackground( Constants.BUTTON_BACKGROUND_COLOR_DEFAULT );
            _click_to_combobox.setForeground( Constants.BUTTON_TEXT_COLOR_DEFAULT );
        }
    }

    void addJButton( final JButton jb, final JPanel p ) {
        jb.setFocusPainted( false );
        jb.setFont( ControlPanel.jcb_font );
        if ( !_configuration.isUseNativeUI() ) {
            jb.setBorder( BorderFactory.createLineBorder( ControlPanel.button_border_color ) );
            jb.setBackground( ControlPanel.button_background_color );
            jb.setForeground( ControlPanel.button_text_color );
        }
        p.add( jb );
        jb.addActionListener( this );
    }

    void addJCheckBox( final JCheckBox jcb, final JPanel p ) {
        jcb.setFocusPainted( false );
        jcb.setFont( ControlPanel.jcb_font );
        if ( !_configuration.isUseNativeUI() ) {
            jcb.setBackground( ControlPanel.jcb_background_color );
            jcb.setForeground( ControlPanel.jcb_text_color );
        }
        p.add( jcb, "Center" );
        jcb.addActionListener( this );
    }

    void addJTextField( final JTextField tf, final JPanel p ) {
        if ( !_configuration.isUseNativeUI() ) {
            tf.setForeground( ControlPanel.background_color );
            tf.setFont( ControlPanel.jcb_font );
        }
        p.add( tf );
        tf.addActionListener( this );
    }

    void deactivateButtonToReturnToSuperTree() {
        _return_to_super_tree.setText( RETURN_TO_SUPER_TREE_TEXT );
        _return_to_super_tree.setForeground( Constants.BUTTON_TEXT_COLOR_DEFAULT );
        _return_to_super_tree.setEnabled( false );
    }

    void displayedPhylogenyMightHaveChanged( final boolean recalc_longest_ext_node_info ) {
        if ( ( _mainpanel != null ) && ( _mainpanel.getCurrentPhylogeny() != null ) ) {
            if ( getOptions().isShowOverview() ) {
                _mainpanel.getCurrentTreePanel().updateOvSizes();
            }
            _mainpanel.getCurrentTreePanel().recalculateMaxDistanceToRoot();
            setVisibilityOfDomainStrucureControls();
            updateDomainStructureEvaluethresholdDisplay();
            _mainpanel.getCurrentTreePanel().calculateScaleDistance();
            _mainpanel.getCurrentTreePanel().calcMaxDepth();
            _mainpanel.adjustJScrollPane();
            if ( recalc_longest_ext_node_info ) {
                _mainpanel.getCurrentTreePanel().initNodeData();
                _mainpanel.getCurrentTreePanel().calculateLongestExtNodeInfo();
            }
            _mainpanel.getCurrentTreePanel().repaint();
            // _mainpanel.getCurrentTreePanel().setUpUrtFactors();
        }
    }

    void endClickToOptions() {
        _click_to_combobox.addActionListener( this );
    }

    /**
     * Indicates what action should be execute when a node is clicked
     * 
     * @return the click-on action
     */
    NodeClickAction getActionWhenNodeClicked() {
        return _action_when_node_clicked;
    }

    Map<Integer, String> getAllClickToItems() {
        return _all_click_to_names;
    }

    Map<String, Color> getAnnotationColors() {
        return _annotation_colors;
    }

    public JCheckBox getColorBranchesCb() {
        return _color_branches_cb;
    }

    Configuration getConfiguration() {
        return _configuration;
    }

    TreePanel getCurrentTreePanel() {
        return getMainPanel().getCurrentTreePanel();
    }

    public JCheckBox getDisplayAsPhylogramCb() {
        return _display_as_phylogram_cb;
    }

    public JCheckBox getDynamicallyHideData() {
        return _dynamically_hide_data;
    }

    private List<Boolean> getIsDrawPhylogramList() {
        return _draw_phylogram;
    }

    MainPanel getMainPanel() {
        return _mainpanel;
    }

    public JCheckBox getNodeDescPopupCb() {
        return _node_desc_popup_cb;
    }

    Options getOptions() {
        return getMainPanel().getOptions();
    }

    private JLabel getSearchFoundCountsLabel() {
        return _search_found_label;
    }

    private JButton getSearchResetButton() {
        return _search_reset_button;
    }

    JTextField getSearchTextField() {
        return _search_tf;
    }

    public JCheckBox getShowEventsCb() {
        return _show_events;
    }

    List<String> getSingleClickToNames() {
        return _click_to_names;
    }

    Map<String, Color> getSpeciesColors() {
        return _species_colors;
    }

    public JCheckBox getWriteConfidenceCb() {
        return _write_confidence;
    }

    private void init() {
        _draw_phylogram = new ArrayList<Boolean>();
        setSpeciesColors( new HashMap<String, Color>() );
        setAnnotationColors( new HashMap<String, Color>() );
    }

    boolean isAntialiasScreenText() {
        return true;
    }

    boolean isColorAccordingToAnnotation() {
        return ( ( _color_according_to_annotation != null ) && _color_according_to_annotation.isSelected() );
    }

    boolean isColorAccordingToTaxonomy() {
        return ( ( _color_acc_species != null ) && _color_acc_species.isSelected() );
    }

    boolean isColorBranches() {
        return ( ( ( getColorBranchesCb() != null ) && getColorBranchesCb().isSelected() ) || ( ( getColorBranchesCb() == null ) && _color_branches ) );
    }

    boolean isDrawPhylogram() {
        return isDrawPhylogram( getMainPanel().getCurrentTabIndex() );
    }

    private boolean isDrawPhylogram( final int index ) {
        return getIsDrawPhylogramList().get( index );
    }

    boolean isDynamicallyHideData() {
        return ( ( getDynamicallyHideData() != null ) && getDynamicallyHideData().isSelected() );
    }

    boolean isEvents() {
        return ( ( getShowEventsCb() != null ) && getShowEventsCb().isSelected() );
    }

    boolean isNodeDescPopup() {
        return ( ( getNodeDescPopupCb() != null ) && getNodeDescPopupCb().isSelected() );
    }

    boolean isShowAnnotation() {
        return ( ( _show_annotation != null ) && _show_annotation.isSelected() );
    }

    boolean isShowBinaryCharacterCounts() {
        return ( ( _show_binary_character_counts != null ) && _show_binary_character_counts.isSelected() );
    }

    boolean isShowBinaryCharacters() {
        return ( ( _show_binary_characters != null ) && _show_binary_characters.isSelected() );
    }

    boolean isShowBootstrapValues() {
        return ( ( getWriteConfidenceCb() != null ) && getWriteConfidenceCb().isSelected() );
    }

    boolean isShowDate() {
        // TODO Auto-generated method stub
        return true;
    }

    boolean isShowDistribution() {
        return true;
    }

    boolean isShowDomainArchitectures() {
        return ( ( _show_domain_architectures != null ) && _show_domain_architectures.isSelected() );
    }

    boolean isShowGeneNames() {
        return ( ( _show_gene_names != null ) && _show_gene_names.isSelected() );
    }

    boolean isShowGeneSymbols() {
        return ( ( _show_gene_symbols != null ) && _show_gene_symbols.isSelected() );
    }

    boolean isShowInternalData() {
        return ( ( _display_internal_data == null ) || _display_internal_data.isSelected() );
    }

    boolean isShowNodeNames() {
        return ( ( _show_node_names != null ) && _show_node_names.isSelected() );
    }

    public boolean isShowProperties() {
        // TODO Auto-generated method stub
        return false;
    }

    boolean isShowProperty() {
        return ( ( _show_annotation != null ) && _show_annotation.isSelected() );
    }

    boolean isShowSequenceAcc() {
        return ( ( _show_sequence_acc != null ) && _show_sequence_acc.isSelected() );
    }

    boolean isShowTaxonomyCode() {
        return ( ( _show_taxo_code != null ) && _show_taxo_code.isSelected() );
    }

    boolean isShowTaxonomyNames() {
        return ( ( _show_taxo_names != null ) && _show_taxo_names.isSelected() );
    }

    boolean isWidthBranches() {
        return ( ( _width_branches != null ) && _width_branches.isSelected() );
    }

    void phylogenyAdded( final Configuration configuration ) {
        getIsDrawPhylogramList().add( configuration.isDrawAsPhylogram() );
    }

    void phylogenyRemoved( final int index ) {
        getIsDrawPhylogramList().remove( index );
    }

    void search() {
        final MainPanel main_panel = getMainPanel();
        final Phylogeny tree = main_panel.getCurrentPhylogeny();
        if ( ( tree == null ) || tree.isEmpty() ) {
            return;
        }
        String query = getSearchTextField().getText();
        if ( query != null ) {
            query = query.trim();
        }
        else {
            getSearchFoundCountsLabel().setVisible( false );
            getSearchResetButton().setEnabled( false );
            getSearchResetButton().setVisible( false );
            searchReset();
        }
        if ( !ForesterUtil.isEmpty( query ) ) {
            search( main_panel, tree, query );
        }
        else {
            getSearchFoundCountsLabel().setVisible( false );
            getSearchResetButton().setEnabled( false );
            getSearchResetButton().setVisible( false );
            searchReset();
        }
    }

    private void search( final MainPanel atv_panel, final Phylogeny tree, final String query_str ) {
        getSearchFoundCountsLabel().setVisible( true );
        getSearchResetButton().setEnabled( true );
        getSearchResetButton().setVisible( true );
        String[] queries = null;
        Set<PhylogenyNode> nodes = null;
        if ( query_str.indexOf( ',' ) >= 0 ) {
            queries = query_str.split( ",+" );
        }
        else {
            queries = new String[ 1 ];
            queries[ 0 ] = query_str.trim();
        }
        if ( ( queries != null ) && ( queries.length > 0 ) ) {
            nodes = new HashSet<PhylogenyNode>();
            for( String query : queries ) {
                if ( ForesterUtil.isEmpty( query ) ) {
                    continue;
                }
                query = query.trim();
                if ( query.indexOf( '+' ) >= 0 ) {
                    nodes.addAll( PhylogenyMethods.searchDataLogicalAnd( query.split( "\\++" ), tree, getOptions()
                            .isSearchCaseSensitive(), !getOptions().isMatchWholeTermsOnly() ) );
                }
                else {
                    nodes.addAll( PhylogenyMethods.searchData( query,
                                                               tree,
                                                               getOptions().isSearchCaseSensitive(),
                                                               !getOptions().isMatchWholeTermsOnly() ) );
                }
            }
            if ( getOptions().isInverseSearchResult() ) {
                final Set<PhylogenyNode> all = PhylogenyMethods.obtainAllNodesAsSet( tree );
                all.removeAll( nodes );
                nodes = all;
            }
        }
        if ( ( nodes != null ) && ( nodes.size() > 0 ) ) {
            atv_panel.getCurrentTreePanel().setFoundNodes( nodes );
            setSearchFoundCountsOnLabel( nodes.size() );
        }
        else {
            setSearchFoundCountsOnLabel( 0 );
            searchReset();
        }
    }

    private void searchReset() {
        if ( getMainPanel().getCurrentTreePanel() != null ) {
            getMainPanel().getCurrentTreePanel().setFoundNodes( null );
        }
    }

    void setActionWhenNodeClicked( final NodeClickAction action ) {
        _action_when_node_clicked = action;
    }

    void setAnnotationColors( final Map<String, Color> annotation_colors ) {
        _annotation_colors = annotation_colors;
    }

    void setCheckbox( final int which, final boolean state ) {
        switch ( which ) {
            case Configuration.display_as_phylogram:
                if ( getDisplayAsPhylogramCb() != null ) {
                    getDisplayAsPhylogramCb().setSelected( state );
                }
                break;
            case Configuration.display_internal_data:
                if ( _display_internal_data != null ) {
                    _display_internal_data.setSelected( state );
                }
                break;
            case Configuration.color_according_to_species:
                if ( _color_acc_species != null ) {
                    _color_acc_species.setSelected( state );
                }
                break;
            case Configuration.color_according_to_annotation:
                if ( _color_according_to_annotation != null ) {
                    _color_according_to_annotation.setSelected( state );
                }
                break;
            case Configuration.show_node_names:
                if ( _show_node_names != null ) {
                    _show_node_names.setSelected( state );
                }
                break;
            case Configuration.show_tax_code:
                if ( _show_taxo_code != null ) {
                    _show_taxo_code.setSelected( state );
                }
                break;
            case Configuration.show_taxonomy_names:
                if ( _show_taxo_names != null ) {
                    _show_taxo_names.setSelected( state );
                }
                break;
            case Configuration.show_annotation:
                if ( _show_annotation != null ) {
                    _show_annotation.setSelected( state );
                }
                break;
            case Configuration.show_binary_characters:
                if ( _show_binary_characters != null ) {
                    _show_binary_characters.setSelected( state );
                }
                break;
            case Configuration.show_binary_character_counts:
                if ( _show_binary_character_counts != null ) {
                    _show_binary_character_counts.setSelected( state );
                }
                break;
            case Configuration.write_confidence_values:
                if ( getWriteConfidenceCb() != null ) {
                    getWriteConfidenceCb().setSelected( state );
                }
                break;
            case Configuration.write_events:
                if ( getShowEventsCb() != null ) {
                    getShowEventsCb().setSelected( state );
                }
                break;
            case Configuration.color_branches:
                if ( getColorBranchesCb() != null ) {
                    getColorBranchesCb().setSelected( state );
                }
                break;
            case Configuration.width_branches:
                if ( _width_branches != null ) {
                    _width_branches.setSelected( state );
                }
                break;
            case Configuration.show_domain_architectures:
                if ( _show_domain_architectures != null ) {
                    _show_domain_architectures.setSelected( state );
                }
                break;
            case Configuration.show_gene_names:
                if ( _show_gene_names != null ) {
                    _show_gene_names.setSelected( state );
                }
                break;
            case Configuration.show_gene_symbols:
                if ( _show_gene_symbols != null ) {
                    _show_gene_symbols.setSelected( state );
                }
                break;
            case Configuration.show_sequence_acc:
                if ( _show_sequence_acc != null ) {
                    _show_sequence_acc.setSelected( state );
                }
                break;
            case Configuration.dynamically_hide_data:
                if ( getDynamicallyHideData() != null ) {
                    getDynamicallyHideData().setSelected( state );
                }
                break;
            case Configuration.node_data_popup:
                if ( getNodeDescPopupCb() != null ) {
                    getNodeDescPopupCb().setSelected( state );
                }
                break;
            default:
                throw new AssertionError( "unknown checkbox: " + which );
        }
    }

    /**
     * Set this checkbox state. Not all checkboxes have been instantiated
     * depending on the config.
     */
    void setCheckbox( final JCheckBox cb, final boolean state ) {
        if ( cb != null ) {
            cb.setSelected( state );
        }
    }

    void setClickToAction( final int action ) {
        // Set click-to action
        if ( action == _show_data_item ) {
            setActionWhenNodeClicked( NodeClickAction.SHOW_DATA );
        }
        else if ( action == _collapse_cb_item ) {
            setActionWhenNodeClicked( NodeClickAction.COLLAPSE );
        }
        else if ( action == _reroot_cb_item ) {
            setActionWhenNodeClicked( NodeClickAction.REROOT );
        }
        else if ( action == _subtree_cb_item ) {
            setActionWhenNodeClicked( NodeClickAction.SUBTREE );
        }
        else if ( action == _swap_cb_item ) {
            setActionWhenNodeClicked( NodeClickAction.SWAP );
        }
        else if ( action == _color_subtree_cb_item ) {
            setActionWhenNodeClicked( NodeClickAction.COLOR_SUBTREE );
        }
        else if ( action == _open_seq_web_item ) {
            setActionWhenNodeClicked( NodeClickAction.OPEN_SEQ_WEB );
        }
        else if ( action == _blast_item ) {
            if ( !Constants.__RELEASE && !Constants.__SNAPSHOT_RELEASE ) {
                setActionWhenNodeClicked( NodeClickAction.BLAST );
            }
        }
        else if ( action == _open_tax_web_item ) {
            setActionWhenNodeClicked( NodeClickAction.OPEN_TAX_WEB );
        }
        else if ( action == _cut_subtree_item ) {
            setActionWhenNodeClicked( NodeClickAction.CUT_SUBTREE );
        }
        else if ( action == _copy_subtree_item ) {
            setActionWhenNodeClicked( NodeClickAction.COPY_SUBTREE );
        }
        else if ( action == _delete_node_or_subtree_item ) {
            setActionWhenNodeClicked( NodeClickAction.DELETE_NODE_OR_SUBTREE );
        }
        else if ( action == _paste_subtree_item ) {
            setActionWhenNodeClicked( NodeClickAction.PASTE_SUBTREE );
        }
        else if ( action == _add_new_node_item ) {
            setActionWhenNodeClicked( NodeClickAction.ADD_NEW_NODE );
        }
        else if ( action == _edit_node_data_item ) {
            setActionWhenNodeClicked( NodeClickAction.EDIT_NODE_DATA );
        }
        else {
            throw new IllegalStateException( "unknown action: " + action );
        }
        // make sure drop down is displaying the correct action
        // in case this was called from outside the class
        _click_to_combobox.setSelectedIndex( action );
    }

    void setColorBranches( final boolean color_branches ) {
        _color_branches = color_branches;
    }

    /**
     * Sets the color to default values.
     */
    void setDefaultColors() {
        ControlPanel.background_color = Constants.GUI_BACKGROUND_DEFAULT;
        ControlPanel.jcb_text_color = Constants.CHECKBOX_TEXT_COLOR_DEFAULT;
        ControlPanel.jcb_background_color = Constants.CHECKBOX_BACKGROUND_COLOR_DEFAULT;
        ControlPanel.button_text_color = Constants.BUTTON_TEXT_COLOR_DEFAULT;
        ControlPanel.button_background_color = Constants.BUTTON_BACKGROUND_COLOR_DEFAULT;
        ControlPanel.button_border_color = Constants.BUTTON_BORDER_COLOR;
    }// setDefaultColors

    void setDrawPhylogram( final boolean b ) {
        getDisplayAsPhylogramCb().setSelected( b );
        setDrawPhylogram( getMainPanel().getCurrentTabIndex(), b );
    }

    private void setDrawPhylogram( final int index, final boolean b ) {
        getIsDrawPhylogramList().set( index, b );
    }

    void setDrawPhylogramEnabled( final boolean b ) {
        getDisplayAsPhylogramCb().setEnabled( b );
    }

    void setDynamicHidingIsOn( final boolean is_on ) {
        //  if ( !_configuration.isUseNativeUI() ) {
        if ( is_on ) {
            getDynamicallyHideData().setForeground( Constants.CHECKBOX_TEXT_COLOR_ON_DEFAULT );
        }
        else {
            if ( !_configuration.isUseNativeUI() ) {
                getDynamicallyHideData().setForeground( Constants.CHECKBOX_TEXT_COLOR_DEFAULT );
            }
            else {
                getDynamicallyHideData().setForeground( Color.BLACK );
            }
        }
        // }
    }

    private void setSearchFoundCountsOnLabel( final int counts ) {
        getSearchFoundCountsLabel().setText( "Found: " + counts );
    }

    void setShowEvents( final boolean show_events ) {
        if ( getShowEventsCb() == null ) {
            _show_events = new JCheckBox( "" );
        }
        getShowEventsCb().setSelected( show_events );
    }

    void setSpeciesColors( final Map<String, Color> species_colors ) {
        _species_colors = species_colors;
    }

    private void setupClickToOptions() {
        final int default_option = _configuration.getDefaultDisplayClicktoOption();
        int selected_index = 0;
        int cb_index = 0;
        if ( _configuration.doDisplayClickToOption( Configuration.display_node_data ) ) {
            _show_data_item = cb_index;
            addClickToOption( Configuration.display_node_data, _configuration
                    .getClickToTitle( Configuration.display_node_data ) );
            if ( default_option == Configuration.display_node_data ) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if ( _configuration.doDisplayClickToOption( Configuration.collapse_uncollapse ) ) {
            _collapse_cb_item = cb_index;
            addClickToOption( Configuration.collapse_uncollapse, _configuration
                    .getClickToTitle( Configuration.collapse_uncollapse ) );
            if ( default_option == Configuration.collapse_uncollapse ) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if ( _configuration.doDisplayClickToOption( Configuration.reroot ) ) {
            _reroot_cb_item = cb_index;
            addClickToOption( Configuration.reroot, _configuration.getClickToTitle( Configuration.reroot ) );
            if ( default_option == Configuration.reroot ) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if ( _configuration.doDisplayClickToOption( Configuration.subtree ) ) {
            _subtree_cb_item = cb_index;
            addClickToOption( Configuration.subtree, _configuration.getClickToTitle( Configuration.subtree ) );
            if ( default_option == Configuration.subtree ) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if ( _configuration.doDisplayClickToOption( Configuration.swap ) ) {
            _swap_cb_item = cb_index;
            addClickToOption( Configuration.swap, _configuration.getClickToTitle( Configuration.swap ) );
            if ( default_option == Configuration.swap ) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if ( _configuration.doDisplayClickToOption( Configuration.color_subtree ) ) {
            _color_subtree_cb_item = cb_index;
            addClickToOption( Configuration.color_subtree, _configuration.getClickToTitle( Configuration.color_subtree ) );
            if ( default_option == Configuration.color_subtree ) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if ( _configuration.doDisplayClickToOption( Configuration.open_seq_web ) ) {
            _open_seq_web_item = cb_index;
            addClickToOption( Configuration.open_seq_web, _configuration.getClickToTitle( Configuration.open_seq_web ) );
            if ( default_option == Configuration.open_seq_web ) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if ( _configuration.doDisplayClickToOption( Configuration.open_tax_web ) ) {
            _open_tax_web_item = cb_index;
            addClickToOption( Configuration.open_tax_web, _configuration.getClickToTitle( Configuration.open_tax_web ) );
            if ( default_option == Configuration.open_tax_web ) {
                selected_index = cb_index;
            }
            cb_index++;
        }
        if ( getOptions().isEditable() ) {
            if ( _configuration.doDisplayClickToOption( Configuration.cut_subtree ) ) {
                _cut_subtree_item = cb_index;
                addClickToOption( Configuration.cut_subtree, _configuration.getClickToTitle( Configuration.cut_subtree ) );
                if ( default_option == Configuration.cut_subtree ) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if ( _configuration.doDisplayClickToOption( Configuration.copy_subtree ) ) {
                _copy_subtree_item = cb_index;
                addClickToOption( Configuration.copy_subtree, _configuration
                        .getClickToTitle( Configuration.copy_subtree ) );
                if ( default_option == Configuration.copy_subtree ) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if ( _configuration.doDisplayClickToOption( Configuration.paste_subtree ) ) {
                _paste_subtree_item = cb_index;
                addClickToOption( Configuration.paste_subtree, _configuration
                        .getClickToTitle( Configuration.paste_subtree ) );
                if ( default_option == Configuration.paste_subtree ) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if ( _configuration.doDisplayClickToOption( Configuration.delete_subtree_or_node ) ) {
                _delete_node_or_subtree_item = cb_index;
                addClickToOption( Configuration.delete_subtree_or_node, _configuration
                        .getClickToTitle( Configuration.delete_subtree_or_node ) );
                if ( default_option == Configuration.delete_subtree_or_node ) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if ( _configuration.doDisplayClickToOption( Configuration.add_new_node ) ) {
                _add_new_node_item = cb_index;
                addClickToOption( Configuration.add_new_node, _configuration
                        .getClickToTitle( Configuration.add_new_node ) );
                if ( default_option == Configuration.add_new_node ) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if ( _configuration.doDisplayClickToOption( Configuration.edit_node_data ) ) {
                _edit_node_data_item = cb_index;
                addClickToOption( Configuration.edit_node_data, _configuration
                        .getClickToTitle( Configuration.edit_node_data ) );
                if ( default_option == Configuration.edit_node_data ) {
                    selected_index = cb_index;
                }
                cb_index++;
            }
            if ( !Constants.__RELEASE && !Constants.__SNAPSHOT_RELEASE ) {
                if ( _configuration.doDisplayClickToOption( Configuration.blast ) ) {
                    _blast_item = cb_index;
                    addClickToOption( Configuration.blast, _configuration.getClickToTitle( Configuration.blast ) );
                    if ( default_option == Configuration.blast ) {
                        selected_index = cb_index;
                    }
                    cb_index++;
                }
            }
        }
        // Set default selection and its action
        _click_to_combobox.setSelectedIndex( selected_index );
        setClickToAction( selected_index );
    }

    /*
     * Set up the controls from the config settings. 11/26/05
     */
    void setupControls() {
        // The tree display options:
        setupDisplayCheckboxes();
        // Click-to options
        startClickToOptions();
        setupClickToOptions();
        endClickToOptions();
        // Zoom and quick edit buttons
        addButtons();
        setupSearchTools();
    }

    void setUpControlsForDomainStrucures() {
        _domain_display_label = new JLabel( "Domain Display:" );
        add( customizeLabel( _domain_display_label, getConfiguration() ) );
        add( _domain_display_label );
        _zoom_in_domain_structure = new JButton( "d+" );
        _zoom_out_domain_structure = new JButton( "d-" );
        _decr_domain_structure_evalue_thr = new JButton( "-" );
        _incr_domain_structure_evalue_thr = new JButton( "+" );
        _zoom_in_domain_structure.setPreferredSize( new Dimension( 10, 10 ) );
        _zoom_out_domain_structure.setPreferredSize( new Dimension( 10, 10 ) );
        _decr_domain_structure_evalue_thr.setPreferredSize( new Dimension( 10, 10 ) );
        _incr_domain_structure_evalue_thr.setPreferredSize( new Dimension( 10, 10 ) );
        _incr_domain_structure_evalue_thr.setToolTipText( "Increase the E-value threshold by a factor of 10" );
        _decr_domain_structure_evalue_thr.setToolTipText( "Decrease the E-value threshold by a factor of 10" );
        _domain_structure_evalue_thr_tf = new JTextField( 3 );
        _domain_structure_evalue_thr_tf.setEditable( false );
        if ( !getConfiguration().isUseNativeUI() ) {
            _domain_structure_evalue_thr_tf.setForeground( Constants.MENU_BACKGROUND_COLOR_DEFAULT );
            _domain_structure_evalue_thr_tf.setBackground( ControlPanel.jcb_text_color );
            _domain_structure_evalue_thr_tf.setBorder( null );
        }
        final JPanel d1_panel = new JPanel( new GridLayout( 1, 2, 0, 0 ) );
        final JPanel d2_panel = new JPanel( new GridLayout( 1, 3, 0, 0 ) );
        if ( !_configuration.isUseNativeUI() ) {
            d1_panel.setBackground( getBackground() );
            d2_panel.setBackground( getBackground() );
        }
        add( d1_panel );
        add( d2_panel );
        addJButton( _zoom_out_domain_structure, d1_panel );
        addJButton( _zoom_in_domain_structure, d1_panel );
        addJButton( _decr_domain_structure_evalue_thr, d2_panel );
        addJTextField( _domain_structure_evalue_thr_tf, d2_panel );
        addJButton( _incr_domain_structure_evalue_thr, d2_panel );
    }

    private void setupDisplayCheckboxes() {
        if ( _configuration.doDisplayOption( Configuration.display_as_phylogram ) ) {
            addCheckbox( Configuration.display_as_phylogram, _configuration
                    .getDisplayTitle( Configuration.display_as_phylogram ) );
            setCheckbox( Configuration.display_as_phylogram, _configuration
                    .doCheckOption( Configuration.display_as_phylogram ) );
        }
        if ( _configuration.doDisplayOption( Configuration.dynamically_hide_data ) ) {
            addCheckbox( Configuration.dynamically_hide_data, _configuration
                    .getDisplayTitle( Configuration.dynamically_hide_data ) );
            setCheckbox( Configuration.dynamically_hide_data, _configuration
                    .doCheckOption( Configuration.dynamically_hide_data ) );
        }
        if ( _configuration.doDisplayOption( Configuration.node_data_popup ) ) {
            addCheckbox( Configuration.node_data_popup, _configuration.getDisplayTitle( Configuration.node_data_popup ) );
            setCheckbox( Configuration.node_data_popup, _configuration.doCheckOption( Configuration.node_data_popup ) );
        }
        if ( _configuration.doDisplayOption( Configuration.display_internal_data ) ) {
            addCheckbox( Configuration.display_internal_data, _configuration
                    .getDisplayTitle( Configuration.display_internal_data ) );
            setCheckbox( Configuration.display_internal_data, _configuration
                    .doCheckOption( Configuration.display_internal_data ) );
        }
        if ( _configuration.doDisplayOption( Configuration.color_according_to_species ) ) {
            addCheckbox( Configuration.color_according_to_species, _configuration
                    .getDisplayTitle( Configuration.color_according_to_species ) );
            setCheckbox( Configuration.color_according_to_species, _configuration
                    .doCheckOption( Configuration.color_according_to_species ) );
        }
        if ( _configuration.doDisplayOption( Configuration.color_according_to_annotation ) ) {
            addCheckbox( Configuration.color_according_to_annotation, _configuration
                    .getDisplayTitle( Configuration.color_according_to_annotation ) );
            setCheckbox( Configuration.color_according_to_annotation, _configuration
                    .doCheckOption( Configuration.color_according_to_annotation ) );
        }
        if ( _configuration.doDisplayOption( Configuration.color_branches ) ) {
            addCheckbox( Configuration.color_branches, _configuration.getDisplayTitle( Configuration.color_branches ) );
            setCheckbox( Configuration.color_branches, _configuration.doCheckOption( Configuration.color_branches ) );
        }
        if ( _configuration.doDisplayOption( Configuration.width_branches ) ) {
            addCheckbox( Configuration.width_branches, _configuration.getDisplayTitle( Configuration.width_branches ) );
            setCheckbox( Configuration.width_branches, _configuration.doCheckOption( Configuration.width_branches ) );
        }
        final JLabel label = new JLabel( "Display Data:" );
        label.setFont( ControlPanel.jcb_bold_font );
        if ( !getConfiguration().isUseNativeUI() ) {
            label.setForeground( ControlPanel.jcb_text_color );
        }
        add( label );
        if ( _configuration.doDisplayOption( Configuration.show_node_names ) ) {
            addCheckbox( Configuration.show_node_names, _configuration.getDisplayTitle( Configuration.show_node_names ) );
            setCheckbox( Configuration.show_node_names, _configuration.doCheckOption( Configuration.show_node_names ) );
        }
        if ( _configuration.doDisplayOption( Configuration.show_tax_code ) ) {
            addCheckbox( Configuration.show_tax_code, _configuration.getDisplayTitle( Configuration.show_tax_code ) );
            setCheckbox( Configuration.show_tax_code, _configuration.doCheckOption( Configuration.show_tax_code ) );
        }
        if ( _configuration.doDisplayOption( Configuration.show_taxonomy_names ) ) {
            addCheckbox( Configuration.show_taxonomy_names, _configuration
                    .getDisplayTitle( Configuration.show_taxonomy_names ) );
            setCheckbox( Configuration.show_taxonomy_names, _configuration
                    .doCheckOption( Configuration.show_taxonomy_names ) );
        }
        if ( _configuration.doDisplayOption( Configuration.show_gene_symbols ) ) {
            addCheckbox( Configuration.show_gene_symbols, _configuration
                    .getDisplayTitle( Configuration.show_gene_symbols ) );
            setCheckbox( Configuration.show_gene_symbols, _configuration
                    .doCheckOption( Configuration.show_gene_symbols ) );
        }
        if ( _configuration.doDisplayOption( Configuration.show_gene_names ) ) {
            addCheckbox( Configuration.show_gene_names, _configuration.getDisplayTitle( Configuration.show_gene_names ) );
            setCheckbox( Configuration.show_gene_names, _configuration.doCheckOption( Configuration.show_gene_names ) );
        }
        if ( _configuration.doDisplayOption( Configuration.show_sequence_acc ) ) {
            addCheckbox( Configuration.show_sequence_acc, _configuration
                    .getDisplayTitle( Configuration.show_sequence_acc ) );
            setCheckbox( Configuration.show_sequence_acc, _configuration
                    .doCheckOption( Configuration.show_sequence_acc ) );
        }
        if ( _configuration.doDisplayOption( Configuration.show_annotation ) ) {
            addCheckbox( Configuration.show_annotation, _configuration.getDisplayTitle( Configuration.show_annotation ) );
            setCheckbox( Configuration.show_annotation, _configuration.doCheckOption( Configuration.show_annotation ) );
        }
        if ( _configuration.doDisplayOption( Configuration.show_binary_characters ) ) {
            addCheckbox( Configuration.show_binary_characters, _configuration
                    .getDisplayTitle( Configuration.show_binary_characters ) );
            setCheckbox( Configuration.show_binary_characters, _configuration
                    .doCheckOption( Configuration.show_binary_characters ) );
        }
        if ( _configuration.doDisplayOption( Configuration.show_binary_character_counts ) ) {
            addCheckbox( Configuration.show_binary_character_counts, _configuration
                    .getDisplayTitle( Configuration.show_binary_character_counts ) );
            setCheckbox( Configuration.show_binary_character_counts, _configuration
                    .doCheckOption( Configuration.show_binary_character_counts ) );
        }
        if ( _configuration.doDisplayOption( Configuration.show_domain_architectures ) ) {
            addCheckbox( Configuration.show_domain_architectures, _configuration
                    .getDisplayTitle( Configuration.show_domain_architectures ) );
            setCheckbox( Configuration.show_domain_architectures, _configuration
                    .doCheckOption( Configuration.show_domain_architectures ) );
        }
        if ( _configuration.doDisplayOption( Configuration.write_confidence_values ) ) {
            addCheckbox( Configuration.write_confidence_values, _configuration
                    .getDisplayTitle( Configuration.write_confidence_values ) );
            setCheckbox( Configuration.write_confidence_values, _configuration
                    .doCheckOption( Configuration.write_confidence_values ) );
        }
        if ( _configuration.doDisplayOption( Configuration.write_events ) ) {
            addCheckbox( Configuration.write_events, _configuration.getDisplayTitle( Configuration.write_events ) );
            setCheckbox( Configuration.write_events, _configuration.doCheckOption( Configuration.write_events ) );
        }
    }

    void setupSearchTools() {
        final String tip = "Enter text to search for. Use ',' for multiple searches (logical OR) and '+' for logical AND.";
        final JLabel search_label = new JLabel( "Search:" );
        search_label.setFont( ControlPanel.jcb_bold_font );
        if ( !getConfiguration().isUseNativeUI() ) {
            search_label.setForeground( ControlPanel.jcb_text_color );
        }
        add( search_label );
        search_label.setToolTipText( tip );
        _search_found_label = new JLabel();
        getSearchFoundCountsLabel().setVisible( false );
        _search_found_label.setFont( ControlPanel.jcb_bold_font );
        if ( !getConfiguration().isUseNativeUI() ) {
            _search_found_label.setForeground( ControlPanel.jcb_text_color );
        }
        _search_tf = new JTextField( 3 );
        _search_tf.setToolTipText( tip );
        _search_tf.setEditable( true );
        if ( !getConfiguration().isUseNativeUI() ) {
            _search_tf.setForeground( Constants.MENU_BACKGROUND_COLOR_DEFAULT );
            _search_tf.setBackground( ControlPanel.jcb_text_color );
            _search_tf.setBorder( null );
        }
        _search_reset_button = new JButton();
        getSearchResetButton().setText( "Reset" );
        getSearchResetButton().setEnabled( false );
        getSearchResetButton().setVisible( false );
        final JPanel s_panel_1 = new JPanel( new BorderLayout() );
        final JPanel s_panel_2 = new JPanel( new GridLayout( 1, 2, 0, 0 ) );
        s_panel_1.setBackground( getBackground() );
        add( s_panel_1 );
        s_panel_2.setBackground( getBackground() );
        add( s_panel_2 );
        final KeyAdapter key_adapter = new KeyAdapter() {

            @Override
            public void keyReleased( final KeyEvent key_event ) {
                search();
                displayedPhylogenyMightHaveChanged( true );
            }
        };
        final ActionListener action_listener = new ActionListener() {

            @Override
            public void actionPerformed( final ActionEvent e ) {
                searchReset();
                setSearchFoundCountsOnLabel( 0 );
                getSearchFoundCountsLabel().setVisible( false );
                getSearchTextField().setText( "" );
                getSearchResetButton().setEnabled( false );
                getSearchResetButton().setVisible( false );
                displayedPhylogenyMightHaveChanged( true );
            }
        };
        _search_reset_button.addActionListener( action_listener );
        _search_tf.addKeyListener( key_adapter );
        addJTextField( _search_tf, s_panel_1 );
        s_panel_2.add( _search_found_label );
        addJButton( _search_reset_button, s_panel_2 );
    }

    private void setVisibilityOfDomainStrucureControls() {
        if ( _zoom_in_domain_structure != null ) {
            if ( isShowDomainArchitectures() ) {
                _domain_display_label.setVisible( true );
                _zoom_in_domain_structure.setVisible( true );
                _zoom_out_domain_structure.setVisible( true );
                _decr_domain_structure_evalue_thr.setVisible( true );
                _incr_domain_structure_evalue_thr.setVisible( true );
                _domain_structure_evalue_thr_tf.setVisible( true );
            }
            else {
                _domain_display_label.setVisible( false );
                _zoom_in_domain_structure.setVisible( false );
                _zoom_out_domain_structure.setVisible( false );
                _decr_domain_structure_evalue_thr.setVisible( false );
                _incr_domain_structure_evalue_thr.setVisible( false );
                _domain_structure_evalue_thr_tf.setVisible( false );
            }
        }
    }

    /**
     * Fit entire tree into window.
     */
    void showWhole() {
        if ( _mainpanel.getCurrentScrollPane() == null ) {
            return;
        }
        displayedPhylogenyMightHaveChanged( false );
        _mainpanel.getCurrentTreePanel().updateOvSettings();
        _mainpanel.getCurrentTreePanel().validate();
        _mainpanel.validate();
        _mainpanel.getCurrentTreePanel().setParametersForPainting( _mainpanel.getSizeOfViewport().width,
                                                                   _mainpanel.getSizeOfViewport().height,
                                                                   true );
        _mainpanel.getCurrentTreePanel().resetPreferredSize();
        _mainpanel.adjustJScrollPane();
        _mainpanel.getCurrentTreePanel().repaint();
        _mainpanel.getCurrentTreePanel().validate();
        _mainpanel.validate();
        _mainpanel.getCurrentTreePanel().setParametersForPainting( _mainpanel.getSizeOfViewport().width,
                                                                   _mainpanel.getSizeOfViewport().height,
                                                                   true );
        _mainpanel.getCurrentTreePanel().resetPreferredSize();
        _mainpanel.adjustJScrollPane();
        _mainpanel.getCurrentTreePanel().repaint();
        _mainpanel.getCurrentTreePanel().updateOvSizes();
    }

    void showWholeAll() {
        for( final TreePanel tree_panel : _mainpanel.getTreePanels() ) {
            if ( tree_panel != null ) {
                tree_panel.validate();
                tree_panel.setParametersForPainting( _mainpanel.getSizeOfViewport().width, _mainpanel
                        .getSizeOfViewport().height, true );
                tree_panel.resetPreferredSize();
                tree_panel.repaint();
            }
        }
    }

    // Create header for click-to combo box.
    void startClickToOptions() {
        final JLabel spacer = new JLabel( "" );
        spacer.setFont( ControlPanel.jcb_font );
        add( spacer );
        _click_to_label = new JLabel( "Click on Node to:" );
        add( customizeLabel( _click_to_label, getConfiguration() ) );
        _click_to_combobox = new JComboBox();
        _click_to_combobox.setFocusable( false );
        _click_to_combobox.setMaximumRowCount( 14 );
        _click_to_combobox.setFont( ControlPanel.js_font );
        if ( !_configuration.isUseNativeUI() ) {
            _click_to_combobox.setBackground( ControlPanel.background_color );
        }
        // don't add listener until all items are set (or each one will trigger
        // an event)
        // click_to_list.addActionListener(this);
        add( _click_to_combobox );
        // Correlates option names to titles
        _all_click_to_names = new HashMap<Integer, String>();
        _click_to_names = new ArrayList<String>();
    }

    void tabChanged() {
        if ( getMainPanel().getTabbedPane().getTabCount() > 0 ) {
            if ( getCurrentTreePanel().isPhyHasBranchLengths()
                    && ( getCurrentTreePanel().getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                setDrawPhylogramEnabled( true );
                setDrawPhylogram( isDrawPhylogram() );
            }
            else {
                setDrawPhylogramEnabled( false );
                setDrawPhylogram( false );
            }
            if ( getMainPanel().getMainFrame() == null ) {
                // Must be "E" applet version.
                final ArchaeopteryxE e = ( ArchaeopteryxE ) ( ( MainPanelApplets ) getMainPanel() ).getApplet();
                e.setSelectedTypeInTypeMenu( e.getCurrentTreePanel().getPhylogenyGraphicsType() );
            }
            else {
                getMainPanel().getMainFrame().setSelectedTypeInTypeMenu( getMainPanel().getCurrentTreePanel()
                        .getPhylogenyGraphicsType() );
            }
            getMainPanel().getCurrentTreePanel().updateSubSuperTreeButton();
            getMainPanel().getControlPanel().search();
        }
    }

    /**
     * Uncollapse all nodes.
     */
    void uncollapseAll( final TreePanel tp ) {
        final Phylogeny t = tp.getPhylogeny();
        if ( ( t != null ) && !t.isEmpty() ) {
            t.setAllNodesToNotCollapse();
            t.recalculateNumberOfExternalDescendants( false );
            showWhole();
        }
    }

    void updateDomainStructureEvaluethresholdDisplay() {
        if ( _domain_structure_evalue_thr_tf != null ) {
            _domain_structure_evalue_thr_tf.setText( "10^"
                    + _mainpanel.getCurrentTreePanel().getDomainStructureEvalueThreshold() );
        }
    }

    void zoomInX( final float factor, final float x_correction_factor ) {
        final JScrollBar sb = _mainpanel.getCurrentScrollPane().getHorizontalScrollBar();
        final TreePanel treepanel = _mainpanel.getCurrentTreePanel();
        treepanel.multiplyUrtFactor( 1f );
        if ( ( treepanel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR )
                || ( treepanel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                || isDrawPhylogram( _mainpanel.getCurrentTabIndex() )
                || ( getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP ) ) {
            final double x = ( sb.getMaximum() - sb.getMinimum() ) / ( sb.getValue() + ( sb.getVisibleAmount() / 2.0 ) );
            treepanel.setXdistance( ( treepanel.getXdistance() * factor ) );
            treepanel.setXcorrectionFactor( ( treepanel.getXcorrectionFactor() * x_correction_factor ) );
            _mainpanel.adjustJScrollPane();
            treepanel.resetPreferredSize();
            _mainpanel.getCurrentScrollPane().getViewport().validate();
            sb.setValue( ForesterUtil.roundToInt( ( ( sb.getMaximum() - sb.getMinimum() ) / x )
                    - ( sb.getVisibleAmount() / 2.0 ) ) );
        }
        else {
            final int x = sb.getMaximum() - sb.getMinimum() - sb.getVisibleAmount() - sb.getValue();
            treepanel.setXdistance( ( treepanel.getXdistance() * factor ) );
            treepanel.setXcorrectionFactor( ( treepanel.getXcorrectionFactor() * x_correction_factor ) );
            _mainpanel.adjustJScrollPane();
            treepanel.resetPreferredSize();
            _mainpanel.getCurrentScrollPane().getViewport().validate();
            sb.setValue( sb.getMaximum() - sb.getMinimum() - x - sb.getVisibleAmount() );
        }
        treepanel.resetPreferredSize();
        treepanel.updateOvSizes();
    }

    void zoomInY( final float factor ) {
        final JScrollBar sb = _mainpanel.getCurrentScrollPane().getVerticalScrollBar();
        final TreePanel treepanel = _mainpanel.getCurrentTreePanel();
        treepanel.multiplyUrtFactor( 1.1f );
        final double x = ( sb.getMaximum() - sb.getMinimum() ) / ( sb.getValue() + ( sb.getVisibleAmount() / 2.0 ) );
        treepanel.setYdistance( ( treepanel.getYdistance() * factor ) );
        _mainpanel.adjustJScrollPane();
        treepanel.resetPreferredSize();
        _mainpanel.getCurrentScrollPane().getViewport().validate();
        sb.setValue( ForesterUtil.roundToInt( ( ( sb.getMaximum() - sb.getMinimum() ) / x )
                - ( sb.getVisibleAmount() / 2.0 ) ) );
        treepanel.resetPreferredSize();
        treepanel.updateOvSizes();
    }

    void zoomOutX( final float factor, final float x_correction_factor ) {
        final TreePanel treepanel = _mainpanel.getCurrentTreePanel();
        treepanel.multiplyUrtFactor( 1f );
        if ( ( treepanel.getXdistance() * factor ) > 0.0 ) {
            final JScrollBar sb = _mainpanel.getCurrentScrollPane().getHorizontalScrollBar();
            if ( ( treepanel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR )
                    || ( treepanel.getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                    || isDrawPhylogram( _mainpanel.getCurrentTabIndex() )
                    || ( getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP ) ) {
                _mainpanel.adjustJScrollPane();
                treepanel.resetPreferredSize();
                _mainpanel.getCurrentScrollPane().getViewport().validate();
                final double x = ( sb.getMaximum() - sb.getMinimum() )
                        / ( sb.getValue() + ( sb.getVisibleAmount() / 2.0 ) );
                treepanel.setXdistance( ( treepanel.getXdistance() * factor ) );
                treepanel.setXcorrectionFactor( ( treepanel.getXcorrectionFactor() * x_correction_factor ) );
                _mainpanel.adjustJScrollPane();
                treepanel.resetPreferredSize();
                _mainpanel.getCurrentScrollPane().getViewport().validate();
                sb.setValue( ForesterUtil.roundToInt( ( ( sb.getMaximum() - sb.getMinimum() ) / x )
                        - ( sb.getVisibleAmount() / 2.0 ) ) );
            }
            else {
                final int x = sb.getMaximum() - sb.getMinimum() - sb.getVisibleAmount() - sb.getValue();
                treepanel.setXdistance( treepanel.getXdistance() * factor );
                treepanel.setXcorrectionFactor( treepanel.getXcorrectionFactor() * x_correction_factor );
                if ( x > 0 ) {
                    _mainpanel.adjustJScrollPane();
                    treepanel.resetPreferredSize();
                    _mainpanel.getCurrentScrollPane().getViewport().validate();
                    sb.setValue( sb.getMaximum() - sb.getMinimum() - x - sb.getVisibleAmount() );
                }
            }
            treepanel.resetPreferredSize();
            treepanel.updateOvSizes();
        }
    }

    void zoomOutY( final float factor ) {
        final TreePanel treepanel = _mainpanel.getCurrentTreePanel();
        treepanel.multiplyUrtFactor( 0.9f );
        if ( ( treepanel.getYdistance() * factor ) > 0.0 ) {
            final JScrollBar sb = _mainpanel.getCurrentScrollPane().getVerticalScrollBar();
            final double x = ( sb.getMaximum() - sb.getMinimum() ) / ( sb.getValue() + ( sb.getVisibleAmount() / 2.0 ) );
            treepanel.setYdistance( ( treepanel.getYdistance() * factor ) );
            _mainpanel.adjustJScrollPane();
            treepanel.resetPreferredSize();
            _mainpanel.getCurrentScrollPane().getViewport().validate();
            sb.setValue( ForesterUtil.roundToInt( ( ( sb.getMaximum() - sb.getMinimum() ) / x )
                    - ( sb.getVisibleAmount() / 2.0 ) ) );
            treepanel.resetPreferredSize();
            treepanel.updateOvSizes();
        }
    }

    static JLabel customizeLabel( final JLabel label, final Configuration configuration ) {
        label.setFont( ControlPanel.jcb_bold_font );
        if ( !configuration.isUseNativeUI() ) {
            label.setForeground( ControlPanel.jcb_text_color );
            label.setBackground( ControlPanel.background_color );
        }
        return label;
    }

    enum NodeClickAction {
        SHOW_DATA,
        COLLAPSE,
        REROOT,
        SUBTREE,
        SWAP,
        COLOR_SUBTREE,
        OPEN_TAX_WEB,
        OPEN_SEQ_WEB,
        CUT_SUBTREE,
        COPY_SUBTREE,
        DELETE_NODE_OR_SUBTREE,
        PASTE_SUBTREE,
        ADD_NEW_NODE,
        EDIT_NODE_DATA,
        BLAST;
    }
}
