// $Id: TreePanel.java,v 1.155 2009/12/30 04:33:44 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.InputEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Arc2D;
import java.awt.geom.CubicCurve2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.QuadCurve2D;
import java.awt.geom.Rectangle2D;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URLEncoder;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.JApplet;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JTextArea;
import javax.swing.Popup;
import javax.swing.PopupFactory;

import org.forester.archaeopteryx.ControlPanel.NodeClickAction;
import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public final class TreePanel extends JPanel implements ActionListener, MouseWheelListener, Printable {

    private static final float              PI                                = ( float ) ( Math.PI );
    private static final double             TWO_PI                            = 2 * Math.PI;
    private static final float              ONEHALF_PI                        = ( float ) ( 1.5 * Math.PI );
    private static final float              HALF_PI                           = ( float ) ( Math.PI / 2.0 );
    private static final float              ANGLE_ROTATION_UNIT               = ( float ) ( Math.PI / 32 );
    private static final short              OV_BORDER                         = 10;
    final static Cursor                     CUT_CURSOR                        = Cursor
                                                                                      .getPredefinedCursor( Cursor.CROSSHAIR_CURSOR );
    final static Cursor                     MOVE_CURSOR                       = Cursor
                                                                                      .getPredefinedCursor( Cursor.MOVE_CURSOR );
    final static Cursor                     ARROW_CURSOR                      = Cursor
                                                                                      .getPredefinedCursor( Cursor.DEFAULT_CURSOR );
    final static Cursor                     HAND_CURSOR                       = Cursor
                                                                                      .getPredefinedCursor( Cursor.HAND_CURSOR );
    final static Cursor                     WAIT_CURSOR                       = Cursor
                                                                                      .getPredefinedCursor( Cursor.WAIT_CURSOR );
    private final static long               serialVersionUID                  = -978349745916505029L;
    private final static int                EURO_D                            = 10;
    private final static String             NODE_POPMENU_NODE_CLIENT_PROPERTY = "node";
    private final static int                MIN_ROOT_LENGTH                   = 3;
    private final static int                BOX_SIZE                          = 4;
    private final static int                HALF_BOX_SIZE                     = TreePanel.BOX_SIZE / 2;
    private final static int                MAX_SUBTREES                      = 100;
    private final static int                MAX_NODE_FRAMES                   = 10;
    private final static int                MOVE                              = 20;
    private final static NumberFormat       FORMATTER_CONFIDENCE;
    private final static NumberFormat       FORMATTER_BRANCH_LENGTH;
    private final static int                WIGGLE                            = 2;
    private final static int                HALF_BOX_SIZE_PLUS_WIGGLE         = HALF_BOX_SIZE + WIGGLE;
    private final static int                LIMIT_FOR_HQ_RENDERING            = 1000;
    // TODO "rendering_hints" was static before. Need to make sure everything is OK with it not
    // being static anymore (02/20/2009).
    private final RenderingHints            _rendering_hints                  = new RenderingHints( RenderingHints.KEY_RENDERING,
                                                                                                    RenderingHints.VALUE_RENDER_DEFAULT );
    private File                            _treefile                         = null;
    private Configuration                   _configuration                    = null;
    private final NodeFrame[]               _node_frames                      = new NodeFrame[ TreePanel.MAX_NODE_FRAMES ];
    private int                             _node_frame_index                 = 0;
    private Phylogeny                       _phylogeny                        = null;
    private final Phylogeny[]               _phylogenies                      = new Phylogeny[ TreePanel.MAX_SUBTREES ];
    private int                             _subtree_index                    = 0;
    private MainPanel                       _main_panel                       = null;
    private Set<PhylogenyNode>              _found_nodes                      = null;
    private PhylogenyNode                   _highlight_node                   = null;
    private JPopupMenu                      _node_popup_menu                  = null;
    private JMenuItem                       _node_popup_menu_items[]          = null;
    private int                             _longest_ext_node_info            = 0;
    private float                           _x_correction_factor              = 0.0f;
    private float                           _ov_x_correction_factor           = 0.0f;
    private float                           _x_distance                       = 0.0f;
    private float                           _y_distance                       = 0.0f;
    private PHYLOGENY_GRAPHICS_TYPE         _graphics_type                    = PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR;
    private double                          _domain_structure_width           = Constants.DOMAIN_STRUCTURE_DEFAULT_WIDTH;
    private int                             _domain_structure_e_value_thr_exp = Constants.DOMAIN_STRUCTURE_E_VALUE_THR_DEFAULT_EXP;
    private float                           _last_drag_point_x                = 0;
    private float                           _last_drag_point_y                = 0;
    private ControlPanel                    _control_panel                    = null;
    private int                             _external_node_index              = 0;
    private final Polygon                   _polygon                          = new Polygon();
    private final StringBuilder             _sb                               = new StringBuilder();
    private JColorChooser                   _color_chooser                    = null;
    private double                          _scale_distance                   = 0.0;
    private String                          _scale_label                      = null;
    private final CubicCurve2D              _cubic_curve                      = new CubicCurve2D.Float();
    private final QuadCurve2D               _quad_curve                       = new QuadCurve2D.Float();
    private final Line2D                    _line                             = new Line2D.Float();
    private final Ellipse2D                 _ellipse                          = new Ellipse2D.Float();
    private final Rectangle2D               _rectangle                        = new Rectangle2D.Float();
    private Options                         _options                          = null;
    private float                           _ov_max_width                     = 0;
    private float                           _ov_max_height                    = 0;
    private int                             _ov_x_position                    = 0;
    private int                             _ov_y_position                    = 0;
    private int                             _ov_y_start                       = 0;
    private float                           _ov_y_distance                    = 0;
    private float                           _ov_x_distance                    = 0;
    private boolean                         _ov_on                            = false;
    private double                          _urt_starting_angle               = ( float ) ( Math.PI / 2 );
    private float                           _urt_factor                       = 1;
    private float                           _urt_factor_ov                    = 1;
    private final boolean                   _phy_has_branch_lengths;
    private final Rectangle2D               _ov_rectangle                     = new Rectangle2D.Float();
    private boolean                         _in_ov_rect                       = false;
    private boolean                         _in_ov                            = false;
    private final Rectangle                 _ov_virtual_rectangle             = new Rectangle();
    final private static double             _180_OVER_PI                      = 180.0 / Math.PI;
    private static final float              ROUNDED_D                         = 8;
    private int                             _circ_max_depth;
    private int                             _circ_num_ext_nodes;
    private PhylogenyNode                   _root;
    final private Arc2D                     _arc                              = new Arc2D.Double();
    final private HashMap<Integer, Double>  _urt_nodeid_angle_map             = new HashMap<Integer, Double>();
    final private HashMap<Integer, Integer> _urt_nodeid_index_map             = new HashMap<Integer, Integer>();
    HashMap<Integer, Short>                 _nodeid_dist_to_leaf              = new HashMap<Integer, Short>();
    private AffineTransform                 _at;
    private double                          _max_distance_to_root             = -1;
    private int                             _dynamic_hiding_factor            = 0;
    private boolean                         _edited                           = false;
    private Popup                           _node_desc_popup;
    private JTextArea                       _rollover_popup;
    //private final short                     _skip_counter                     = 0;
    private final StringBuffer              _popup_buffer                     = new StringBuffer();
    final private static Font               POPUP_FONT                        = new Font( Configuration
                                                                                                  .getDefaultFontFamilyName(),
                                                                                          Font.PLAIN,
                                                                                          12 );
    static {
        final DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator( '.' );
        FORMATTER_CONFIDENCE = new DecimalFormat( "#.###", dfs );
        FORMATTER_BRANCH_LENGTH = new DecimalFormat( "#.###", dfs );
    }

    TreePanel( final Phylogeny t, final Configuration configuration, final MainPanel tjp ) {
        requestFocusInWindow();
        addKeyListener( new KeyAdapter() {

            @Override
            public void keyPressed( final KeyEvent key_event ) {
                keyPressedCalls( key_event );
                requestFocusInWindow();
            }
        } );
        addFocusListener( new FocusAdapter() {

            @Override
            public void focusGained( final FocusEvent e ) {
                requestFocusInWindow();
            }
        } );
        if ( ( t == null ) || t.isEmpty() ) {
            throw new IllegalArgumentException( "ill advised attempt to draw phylogeny which is null or empty" );
        }
        _graphics_type = tjp.getOptions().getPhylogenyGraphicsType();
        _main_panel = tjp;
        _configuration = configuration;
        _phylogeny = t;
        _phy_has_branch_lengths = ForesterUtil.isHasAtLeastOneBranchLengthLargerThanZero( _phylogeny );
        init();
        if ( !_phylogeny.isEmpty() ) {
            _phylogeny.recalculateNumberOfExternalDescendants( true );
        }
        setBackground( getTreeColorSet().getBackgroundColor() );
        final MouseListener mouse_listener = new MouseListener( this );
        addMouseListener( mouse_listener );
        addMouseMotionListener( mouse_listener );
        addMouseWheelListener( this );
        calculateScaleDistance();
        FORMATTER_CONFIDENCE.setMaximumFractionDigits( configuration.getNumberOfDigitsAfterCommaForConfidenceValues() );
        FORMATTER_BRANCH_LENGTH.setMaximumFractionDigits( configuration
                .getNumberOfDigitsAfterCommaForBranchLengthValues() );
    }

    final public void actionPerformed( final ActionEvent e ) {
        int index;
        boolean done = false;
        final JMenuItem node_popup_menu_item = ( JMenuItem ) e.getSource();
        for( index = 0; ( index < _node_popup_menu_items.length ) && !done; index++ ) {
            // NOTE: index corresponds to the indices of click-to options
            // in the control panel.
            if ( node_popup_menu_item == _node_popup_menu_items[ index ] ) {
                // Set this as the new default click-to action
                _main_panel.getControlPanel().setClickToAction( index );
                final PhylogenyNode node = ( PhylogenyNode ) _node_popup_menu
                        .getClientProperty( NODE_POPMENU_NODE_CLIENT_PROPERTY );
                handleClickToAction( _control_panel.getActionWhenNodeClicked(), node );
                done = true;
            }
        }
        repaint();
        requestFocusInWindow();
    }

    final private void addEmptyNode( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        final String label = getASimpleTextRepresentationOfANode( node );
        String msg = "";
        if ( ForesterUtil.isEmpty( label ) ) {
            msg = "How to add the new, empty node?";
        }
        else {
            msg = "How to add the new, empty node to node" + label + "?";
        }
        final Object[] options = { "As sibling", "As descendant", "Cancel" };
        final int r = JOptionPane.showOptionDialog( this,
                                                    msg,
                                                    "Addition of Empty New Node",
                                                    JOptionPane.CLOSED_OPTION,
                                                    JOptionPane.QUESTION_MESSAGE,
                                                    null,
                                                    options,
                                                    options[ 2 ] );
        boolean add_as_sibling = true;
        if ( r == 1 ) {
            add_as_sibling = false;
        }
        else if ( r != 0 ) {
            return;
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( new PhylogenyNode() );
        phy.setRooted( true );
        if ( add_as_sibling ) {
            if ( node.isRoot() ) {
                JOptionPane.showMessageDialog( this,
                                               "Cannot add sibling to root",
                                               "Attempt to add sibling to root",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            phy.addAsSibling( node );
        }
        else {
            phy.addAsChild( node );
        }
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.hashIDs();
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        setEdited( true );
        repaint();
    }

    final private void assignGraphicsForBranchWithColorForParentBranch( final PhylogenyNode node,
                                                                        final boolean is_vertical,
                                                                        final Graphics g,
                                                                        final boolean to_pdf,
                                                                        final boolean to_graphics_file ) {
        final NodeClickAction action = _control_panel.getActionWhenNodeClicked();
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else if ( ( ( action == NodeClickAction.COPY_SUBTREE ) || ( action == NodeClickAction.CUT_SUBTREE )
                || ( action == NodeClickAction.DELETE_NODE_OR_SUBTREE ) || ( action == NodeClickAction.PASTE_SUBTREE ) || ( action == NodeClickAction.ADD_NEW_NODE ) )
                && ( getCutOrCopiedTree() != null )
                && ( getCopiedAndPastedNodes() != null )
                && !to_pdf
                && !to_graphics_file && getCopiedAndPastedNodes().contains( node ) ) {
            g.setColor( getTreeColorSet().getFoundColor() );
        }
        else if ( getControlPanel().isColorBranches() && ( PhylogenyMethods.getBranchColorValue( node ) != null ) ) {
            g.setColor( PhylogenyMethods.getBranchColorValue( node ) );
        }
        else if ( to_pdf ) {
            g.setColor( getTreeColorSet().getBranchColorForPdf() );
        }
        else {
            g.setColor( getTreeColorSet().getBranchColor() );
        }
    }

    final void assignGraphicsForNodeBoxWithColorForParentBranch( final PhylogenyNode node, final Graphics g ) {
        if ( getControlPanel().isColorBranches() && ( PhylogenyMethods.getBranchColorValue( node ) != null ) ) {
            g.setColor( PhylogenyMethods.getBranchColorValue( node ) );
        }
        else {
            g.setColor( getTreeColorSet().getBranchColor() );
        }
    }

    final private void blast( final PhylogenyNode node ) {
        if ( !isCanBlast( node ) ) {
            JOptionPane.showMessageDialog( this,
                                           "No sequence information present",
                                           "Cannot Blast",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( node.getNodeData().isHasSequence() ) {
            String name = "";
            if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getName() ) ) {
                name = node.getNodeData().getSequence().getName();
            }
            else if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getSymbol() ) ) {
                name = node.getNodeData().getSequence().getSymbol();
            }
            else if ( node.getNodeData().getSequence().getAccession() != null ) {
                name = node.getNodeData().getSequence().getAccession().getValue();
            }
            if ( !ForesterUtil.isEmpty( name ) ) {
                try {
                    System.out.println( "trying: " + name );
                    final Blast s = new Blast();
                    s.go( name );
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                }
            }
        }
    }

    final void calcMaxDepth() {
        if ( _phylogeny != null ) {
            _circ_max_depth = PhylogenyMethods.calculateMaxDepth( _phylogeny );
        }
    }

    /**
     * Calculate the length of the distance between the given node and its
     * parent.
     * 
     * @param node
     * @param ext_node_x
     * @factor
     * @return the distance value
     */
    final private float calculateBranchLengthToParent( final PhylogenyNode node, final float factor ) {
        if ( getControlPanel().isDrawPhylogram() ) {
            if ( node.getDistanceToParent() < 0.0 ) {
                return 0.0f;
            }
            return ( float ) ( getXcorrectionFactor() * node.getDistanceToParent() );
        }
        else {
            if ( ( factor == 0 ) || isNonLinedUpCladogram() ) {
                return getXdistance();
            }
            return getXdistance() * factor;
        }
    }

    final private Color calculateColorForAnnotation( final PhylogenyData ann ) {
        Color c = getTreeColorSet().getAnnotationColor();
        if ( getControlPanel().isColorAccordingToAnnotation() && ( getControlPanel().getAnnotationColors() != null ) ) {
            c = getControlPanel().getAnnotationColors().get( ann.asSimpleText().toString() );
            if ( c == null ) {
                c = getTreeColorSet().getAnnotationColor();
            }
        }
        return c;
    }

    final void calculateLongestExtNodeInfo() {
        if ( ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return;
        }
        int longest = 20;
        for( final PhylogenyNode node : _phylogeny.getExternalNodes() ) {
            int sum = 0;
            if ( node.isCollapse() ) {
                continue;
            }
            if ( getControlPanel().isShowNodeNames() ) {
                sum += getTreeFontSet()._fm_large.stringWidth( node.getNodeName() + " " );
            }
            if ( node.getNodeData().isHasSequence() ) {
                if ( getControlPanel().isShowSequenceAcc()
                        && ( node.getNodeData().getSequence().getAccession() != null ) ) {
                    sum += getTreeFontSet()._fm_large.stringWidth( node.getNodeData().getSequence().getAccession()
                            .getValue()
                            + " " );
                }
                if ( getControlPanel().isShowGeneNames() && ( node.getNodeData().getSequence().getName().length() > 0 ) ) {
                    sum += getTreeFontSet()._fm_large.stringWidth( node.getNodeData().getSequence().getName() + " " );
                }
                if ( getControlPanel().isShowGeneSymbols()
                        && ( node.getNodeData().getSequence().getSymbol().length() > 0 ) ) {
                    sum += getTreeFontSet()._fm_large.stringWidth( node.getNodeData().getSequence().getSymbol() + " " );
                }
                if ( getControlPanel().isShowAnnotation()
                        && ( node.getNodeData().getSequence().getAnnotations() != null )
                        && !node.getNodeData().getSequence().getAnnotations().isEmpty() ) {
                    sum += getTreeFontSet()._fm_large.stringWidth( node.getNodeData().getSequence().getAnnotations()
                            .get( 0 ).asSimpleText()
                            + " " );
                }
            }
            if ( node.getNodeData().isHasTaxonomy() ) {
                final Taxonomy tax = node.getNodeData().getTaxonomy();
                if ( getControlPanel().isShowTaxonomyCode() && !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                    sum += getTreeFontSet()._fm_large_italic.stringWidth( tax.getTaxonomyCode() + " " );
                }
                if ( getControlPanel().isShowTaxonomyNames() && !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                    sum += getTreeFontSet()._fm_large_italic.stringWidth( tax.getScientificName() + " " );
                }
                if ( getControlPanel().isShowTaxonomyNames() && !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
                    sum += getTreeFontSet()._fm_large_italic.stringWidth( tax.getCommonName() + " " );
                }
            }
            if ( getControlPanel().isShowBinaryCharacters() && node.getNodeData().isHasBinaryCharacters() ) {
                sum += getTreeFontSet()._fm_large.stringWidth( node.getNodeData().getBinaryCharacters()
                        .getGainedCharactersAsStringBuffer().toString() );
            }
            if ( getControlPanel().isShowDomainArchitectures() && node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                sum += ( ( RenderableDomainArchitecture ) node.getNodeData().getSequence().getDomainArchitecture() )
                        .getRenderingSize().getWidth();
            }
            if ( sum >= Constants.EXT_NODE_INFO_LENGTH_MAX ) {
                setLongestExtNodeInfo( Constants.EXT_NODE_INFO_LENGTH_MAX );
                return;
            }
            if ( sum > longest ) {
                longest = sum;
            }
        }
        if ( longest >= Constants.EXT_NODE_INFO_LENGTH_MAX ) {
            setLongestExtNodeInfo( Constants.EXT_NODE_INFO_LENGTH_MAX );
        }
        else {
            setLongestExtNodeInfo( longest );
        }
    }

    final private float calculateOvBranchLengthToParent( final PhylogenyNode node, final int factor ) {
        if ( getControlPanel().isDrawPhylogram() ) {
            if ( node.getDistanceToParent() < 0.0 ) {
                return 0.0f;
            }
            return ( float ) ( getOvXcorrectionFactor() * node.getDistanceToParent() );
        }
        else {
            if ( ( factor == 0 ) || isNonLinedUpCladogram() ) {
                return getOvXDistance();
            }
            return getOvXDistance() * factor;
        }
    }

    final void calculateScaleDistance() {
        if ( ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return;
        }
        final double height = getMaxDistanceToRoot();
        if ( height > 0 ) {
            if ( ( height <= 0.5 ) ) {
                setScaleDistance( 0.01 );
            }
            else if ( height <= 5.0 ) {
                setScaleDistance( 0.1 );
            }
            else if ( height <= 50.0 ) {
                setScaleDistance( 1 );
            }
            else if ( height <= 500.0 ) {
                setScaleDistance( 10 );
            }
            else {
                setScaleDistance( 100 );
            }
        }
        else {
            setScaleDistance( 0.0 );
        }
        String scale_label = String.valueOf( getScaleDistance() );
        if ( !ForesterUtil.isEmpty( _phylogeny.getDistanceUnit() ) ) {
            scale_label += " [" + _phylogeny.getDistanceUnit() + "]";
        }
        setScaleLabel( scale_label );
    }

    final private void cannotOpenBrowserWarningMessage( final String type_type ) {
        JOptionPane.showMessageDialog( this,
                                       "Cannot launch web browser for " + type_type + " data of this node",
                                       "Cannot launch web browser",
                                       JOptionPane.WARNING_MESSAGE );
    }

    /**
     * Collapse the tree from the given node
     * 
     * @param node
     *            a PhylogenyNode
     */
    final void collapse( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot collapse in unrooted display type",
                                           "Attempt to collapse in unrooted display",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( !node.isExternal() && !node.isRoot() ) {
            final boolean collapse = !node.isCollapse();
            Util.collapseSubtree( node, collapse );
            _phylogeny.recalculateNumberOfExternalDescendants( true );
            resetNodeIdToDistToLeafMap();
            calculateLongestExtNodeInfo();
            resetPreferredSize();
            updateOvSizes();
            _main_panel.adjustJScrollPane();
            repaint();
        }
    }

    final void collapseSpeciesSpecificSubtrees() {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        setWaitCursor();
        Util.collapseSpeciesSpecificSubtrees( _phylogeny );
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        calculateLongestExtNodeInfo();
        resetPreferredSize();
        _main_panel.adjustJScrollPane();
        setArrowCursor();
        repaint();
    }

    final private void colorizeSubtree( final Color c, final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot colorize subtree in unrooted display type",
                                           "Attempt to colorize subtree in unrooted display",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        _control_panel.setColorBranches( true );
        if ( _control_panel.getColorBranchesCb() != null ) {
            _control_panel.getColorBranchesCb().setSelected( true );
        }
        for( final PreorderTreeIterator it = new PreorderTreeIterator( node ); it.hasNext(); ) {
            it.next().getBranchData().setBranchColor( new BranchColor( c ) );
        }
        repaint();
    }

    final private void colorSubtree( final PhylogenyNode node ) {
        Color intitial_color = null;
        if ( getControlPanel().isColorBranches() && ( PhylogenyMethods.getBranchColorValue( node ) != null )
                && ( ( ( !node.isRoot() && ( node.getParent().getNumberOfDescendants() < 3 ) ) ) || ( node.isRoot() ) ) ) {
            intitial_color = PhylogenyMethods.getBranchColorValue( node );
        }
        else {
            intitial_color = getTreeColorSet().getBranchColor();
        }
        _color_chooser.setColor( intitial_color );
        _color_chooser.setPreviewPanel( new JPanel() );
        final JDialog dialog = JColorChooser
                .createDialog( this,
                               "Subtree colorization",
                               true,
                               _color_chooser,
                               new SubtreeColorizationActionListener( _color_chooser, node ),
                               null );
        dialog.setVisible( true );
    }

    final void confColor() {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        setWaitCursor();
        Util.colorPhylogenyAccordingToConfidenceValues( _phylogeny, this );
        _control_panel.setColorBranches( true );
        if ( _control_panel.getColorBranchesCb() != null ) {
            _control_panel.getColorBranchesCb().setSelected( true );
        }
        setArrowCursor();
        repaint();
    }

    final private void copySubtree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        setCutOrCopiedTree( _phylogeny.subTree( node ) );
        final Set<PhylogenyNode> nodes = PhylogenyMethods.getAllDescendants( node );
        nodes.add( node );
        setCopiedAndPastedNodes( nodes );
        repaint();
    }

    final private void cutSubtree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        if ( node.isRoot() ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot cut entire tree as subtree",
                                           "Attempt to cut entire tree",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        final String label = getASimpleTextRepresentationOfANode( node );
        final int r = JOptionPane.showConfirmDialog( null,
                                                     "Cut subtree" + label + "?",
                                                     "Confirm Cutting of Subtree",
                                                     JOptionPane.YES_NO_OPTION );
        if ( r != JOptionPane.OK_OPTION ) {
            return;
        }
        setCopiedAndPastedNodes( null );
        setCutOrCopiedTree( _phylogeny.subTree( node ) );
        _phylogeny.deleteSubtree( node, true );
        _phylogeny.hashIDs();
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        setEdited( true );
        repaint();
    }

    final private void cycleColors() {
        getMainPanel().getTreeColorSet().cycleColorScheme();
        for( final TreePanel tree_panel : getMainPanel().getTreePanels() ) {
            tree_panel.setBackground( getMainPanel().getTreeColorSet().getBackgroundColor() );
        }
    }

    final void decreaseDomainStructureEvalueThreshold() {
        if ( _domain_structure_e_value_thr_exp > -20 ) {
            _domain_structure_e_value_thr_exp -= 1;
        }
    }

    final private void decreaseOvSize() {
        if ( ( getOvMaxWidth() > 20 ) && ( getOvMaxHeight() > 20 ) ) {
            setOvMaxWidth( getOvMaxWidth() - 5 );
            setOvMaxHeight( getOvMaxHeight() - 5 );
            updateOvSettings();
            getControlPanel().displayedPhylogenyMightHaveChanged( false );
        }
    }

    final private void deleteNodeOrSubtree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        if ( node.isRoot() ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot delete entire tree",
                                           "Attempt to delete entire tree",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        final String label = getASimpleTextRepresentationOfANode( node );
        final Object[] options = { "Node only", "Entire subtree", "Cancel" };
        final int r = JOptionPane.showOptionDialog( this,
                                                    "Delete" + label + "?",
                                                    "Delete Node/Subtree",
                                                    JOptionPane.CLOSED_OPTION,
                                                    JOptionPane.QUESTION_MESSAGE,
                                                    null,
                                                    options,
                                                    options[ 2 ] );
        boolean node_only = true;
        if ( r == 1 ) {
            node_only = false;
        }
        else if ( r != 0 ) {
            return;
        }
        if ( node_only ) {
            PhylogenyMethods.removeNode( node, _phylogeny );
        }
        else {
            _phylogeny.deleteSubtree( node, true );
        }
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.hashIDs();
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        setEdited( true );
        repaint();
    }

    final private void displayNodePopupMenu( final PhylogenyNode node, final int x, final int y ) {
        makePopupMenus( node );
        _node_popup_menu.putClientProperty( NODE_POPMENU_NODE_CLIENT_PROPERTY, node );
        _node_popup_menu.show( this, x, y );
    }

    final private void drawArc( final double x,
                                final double y,
                                final double width,
                                final double heigth,
                                final double start_angle,
                                final double arc_angle,
                                final Graphics2D g ) {
        _arc.setArc( x, y, width, heigth, _180_OVER_PI * start_angle, _180_OVER_PI * arc_angle, Arc2D.OPEN );
        g.draw( _arc );
    }

    final private void drawLine( final double x1, final double y1, final double x2, final double y2, final Graphics2D g ) {
        if ( ( x1 == x2 ) && ( y1 == y2 ) ) {
            return;
        }
        _line.setLine( x1, y1, x2, y2 );
        g.draw( _line );
    }

    final private void drawOval( final double x,
                                 final double y,
                                 final double width,
                                 final double heigth,
                                 final Graphics2D g ) {
        _ellipse.setFrame( x, y, width, heigth );
        g.draw( _ellipse );
    }

    final private void drawOvalFilled( final double x,
                                       final double y,
                                       final double width,
                                       final double heigth,
                                       final Graphics2D g ) {
        _ellipse.setFrame( x, y, width, heigth );
        g.fill( _ellipse );
    }

    final private void drawRect( final float x, final float y, final float width, final float heigth, final Graphics2D g ) {
        _rectangle.setFrame( x, y, width, heigth );
        g.draw( _rectangle );
    }

    final private void drawRectFilled( final double x,
                                       final double y,
                                       final double width,
                                       final double heigth,
                                       final Graphics2D g ) {
        _rectangle.setFrame( x, y, width, heigth );
        g.fill( _rectangle );
    }

    final private void errorMessageNoCutCopyPasteInUnrootedDisplay() {
        JOptionPane.showMessageDialog( this,
                                       "Cannot cut, copy, paste, add, or delete subtrees/nodes in unrooted display",
                                       "Attempt to cut/copy/paste/add/delete in unrooted display",
                                       JOptionPane.ERROR_MESSAGE );
    }

    /**
     * Find the node, if any, at the given location
     * 
     * @param x
     * @param y
     * @return pointer to the node at x,y, null if not found
     */
    final PhylogenyNode findNode( final int x, final int y ) {
        if ( ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return null;
        }
        for( final PhylogenyNodeIterator iter = _phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( ( _phylogeny.isRooted() || !node.isRoot() || ( node.getNumberOfDescendants() > 2 ) )
                    && ( ( node.getXcoord() - HALF_BOX_SIZE_PLUS_WIGGLE ) <= x )
                    && ( ( node.getXcoord() + HALF_BOX_SIZE_PLUS_WIGGLE ) >= x )
                    && ( ( node.getYcoord() - HALF_BOX_SIZE_PLUS_WIGGLE ) <= y )
                    && ( ( node.getYcoord() + HALF_BOX_SIZE_PLUS_WIGGLE ) >= y ) ) {
                return node;
            }
        }
        return null;
    }

    final private String getASimpleTextRepresentationOfANode( final PhylogenyNode node ) {
        final String tax = PhylogenyMethods.getSpecies( node );
        String label = node.getNodeName();
        if ( !ForesterUtil.isEmpty( label ) && !ForesterUtil.isEmpty( tax ) ) {
            label = label + " " + tax;
        }
        else if ( !ForesterUtil.isEmpty( tax ) ) {
            label = tax;
        }
        else {
            label = "";
        }
        if ( !ForesterUtil.isEmpty( label ) ) {
            label = " [" + label + "]";
        }
        return label;
    }

    final Configuration getConfiguration() {
        return _configuration;
    }

    final ControlPanel getControlPanel() {
        return _control_panel;
    }

    final private Set<PhylogenyNode> getCopiedAndPastedNodes() {
        return getMainPanel().getCopiedAndPastedNodes();
    }

    final private Phylogeny getCutOrCopiedTree() {
        return getMainPanel().getCutOrCopiedTree();
    }

    final int getDomainStructureEvalueThreshold() {
        return _domain_structure_e_value_thr_exp;
    }

    final Set<PhylogenyNode> getFoundNodes() {
        return _found_nodes;
    }

    final private float getLastDragPointX() {
        return _last_drag_point_x;
    }

    final private float getLastDragPointY() {
        return _last_drag_point_y;
    }

    final int getLongestExtNodeInfo() {
        return _longest_ext_node_info;
    }

    final public MainPanel getMainPanel() {
        return _main_panel;
    }

    final private short getMaxBranchesToLeaf( final PhylogenyNode node ) {
        if ( !_nodeid_dist_to_leaf.containsKey( node.getNodeId() ) ) {
            final short m = PhylogenyMethods.calculateMaxBranchesToLeaf( node );
            _nodeid_dist_to_leaf.put( node.getNodeId(), m );
            return m;
        }
        else {
            return _nodeid_dist_to_leaf.get( node.getNodeId() );
        }
    }

    final private double getMaxDistanceToRoot() {
        if ( _max_distance_to_root < 0 ) {
            recalculateMaxDistanceToRoot();
        }
        return _max_distance_to_root;
    }

    final Options getOptions() {
        if ( _options == null ) {
            _options = getControlPanel().getOptions();
        }
        return _options;
    }

    final private float getOvMaxHeight() {
        return _ov_max_height;
    }

    final private float getOvMaxWidth() {
        return _ov_max_width;
    }

    final Rectangle2D getOvRectangle() {
        return _ov_rectangle;
    }

    final Rectangle getOvVirtualRectangle() {
        return _ov_virtual_rectangle;
    }

    final private float getOvXcorrectionFactor() {
        return _ov_x_correction_factor;
    }

    final private float getOvXDistance() {
        return _ov_x_distance;
    }

    final private int getOvXPosition() {
        return _ov_x_position;
    }

    final private float getOvYDistance() {
        return _ov_y_distance;
    }

    final private int getOvYPosition() {
        return _ov_y_position;
    }

    final private int getOvYStart() {
        return _ov_y_start;
    }

    /**
     * Get a pointer to the phylogeny 
     * 
     * @return a pointer to the phylogeny
     */
    final Phylogeny getPhylogeny() {
        return _phylogeny;
    }

    final PHYLOGENY_GRAPHICS_TYPE getPhylogenyGraphicsType() {
        return _graphics_type;
    }

    final private double getScaleDistance() {
        return _scale_distance;
    }

    final private String getScaleLabel() {
        return _scale_label;
    }

    final double getStartingAngle() {
        return _urt_starting_angle;
    }

    /**
     * @return pointer to colorset for tree drawing
     */
    final TreeColorSet getTreeColorSet() {
        return getMainPanel().getTreeColorSet();
    }

    final File getTreeFile() {
        return _treefile;
    }

    final private TreeFontSet getTreeFontSet() {
        return getMainPanel().getTreeFontSet();
    }

    final private float getUrtFactor() {
        return _urt_factor;
    }

    final private float getUrtFactorOv() {
        return _urt_factor_ov;
    }

    final float getXcorrectionFactor() {
        return _x_correction_factor;
    }

    final float getXdistance() {
        return _x_distance;
    }

    final float getYdistance() {
        return _y_distance;
    }

    final private void handleClickToAction( final NodeClickAction action, final PhylogenyNode node ) {
        switch ( action ) {
            case SHOW_DATA:
                showNodeFrame( node );
                break;
            case COLLAPSE:
                collapse( node );
                break;
            case REROOT:
                reRoot( node );
                break;
            case SUBTREE:
                subTree( node );
                break;
            case SWAP:
                swap( node );
                break;
            case COLOR_SUBTREE:
                colorSubtree( node );
                break;
            case OPEN_SEQ_WEB:
                openSeqWeb( node );
                break;
            case BLAST:
                blast( node );
                break;
            case OPEN_TAX_WEB:
                openTaxWeb( node );
                break;
            case CUT_SUBTREE:
                cutSubtree( node );
                break;
            case COPY_SUBTREE:
                copySubtree( node );
                break;
            case PASTE_SUBTREE:
                pasteSubtree( node );
                break;
            case DELETE_NODE_OR_SUBTREE:
                deleteNodeOrSubtree( node );
                break;
            case ADD_NEW_NODE:
                addEmptyNode( node );
                break;
            case EDIT_NODE_DATA:
                showNodeEditFrame( node );
                break;
            default:
                throw new IllegalArgumentException( "unknown action: " + action );
        }
    }

    final void increaseDomainStructureEvalueThreshold() {
        if ( _domain_structure_e_value_thr_exp < 3 ) {
            _domain_structure_e_value_thr_exp += 1;
        }
    }

    final private void increaseOvSize() {
        if ( ( getOvMaxWidth() < getMainPanel().getCurrentScrollPane().getViewport().getVisibleRect().getWidth() / 2 )
                && ( getOvMaxHeight() < getMainPanel().getCurrentScrollPane().getViewport().getVisibleRect()
                        .getHeight() / 2 ) ) {
            setOvMaxWidth( getOvMaxWidth() + 5 );
            setOvMaxHeight( getOvMaxHeight() + 5 );
            updateOvSettings();
            getControlPanel().displayedPhylogenyMightHaveChanged( false );
        }
    }

    final void inferCommonPartOfScientificNames() {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        setWaitCursor();
        Util.inferCommonPartOfScientificNames( _phylogeny );
        setArrowCursor();
        repaint();
    }

    final private void init() {
        _color_chooser = new JColorChooser();
        _rollover_popup = new JTextArea();
        _rollover_popup.setFont( POPUP_FONT );
        resetNodeIdToDistToLeafMap();
        setTextAntialias();
        setTreeFile( null );
        setEdited( false );
        initializeOvSettings();
        setStartingAngle( TWO_PI * 3 / 4 );
    }

    final private void initializeOvSettings() {
        setOvMaxHeight( getConfiguration().getOvMaxHeight() );
        setOvMaxWidth( getConfiguration().getOvMaxWidth() );
    }

    final void initNodeData() {
        if ( ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return;
        }
        double max_original_domain_structure_width = 0.0;
        for( final PhylogenyNode node : _phylogeny.getExternalNodes() ) {
            if ( node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                RenderableDomainArchitecture rds = null;
                if ( !( node.getNodeData().getSequence().getDomainArchitecture() instanceof RenderableDomainArchitecture ) ) {
                    rds = new RenderableDomainArchitecture( node.getNodeData().getSequence().getDomainArchitecture() );
                    node.getNodeData().getSequence().setDomainArchitecture( rds );
                }
                else {
                    rds = ( RenderableDomainArchitecture ) node.getNodeData().getSequence().getDomainArchitecture();
                }
                if ( getControlPanel().isShowDomainArchitectures() ) {
                    final double dsw = rds.getOriginalSize().getWidth();
                    if ( dsw > max_original_domain_structure_width ) {
                        max_original_domain_structure_width = dsw;
                    }
                }
            }
        }
        if ( getControlPanel().isShowDomainArchitectures() ) {
            final double ds_factor_width = _domain_structure_width / max_original_domain_structure_width;
            for( final PhylogenyNode node : _phylogeny.getExternalNodes() ) {
                if ( node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                    final RenderableDomainArchitecture rds = ( RenderableDomainArchitecture ) node.getNodeData()
                            .getSequence().getDomainArchitecture();
                    rds.setRenderingFactorWidth( ds_factor_width );
                    rds.setParameter( _domain_structure_e_value_thr_exp );
                }
            }
        }
    }

    final boolean inOv( final MouseEvent e ) {
        return ( ( e.getX() > getVisibleRect().x + getOvXPosition() + 1 )
                && ( e.getX() < getVisibleRect().x + getOvXPosition() + getOvMaxWidth() - 1 )
                && ( e.getY() > getVisibleRect().y + getOvYPosition() + 1 ) && ( e.getY() < getVisibleRect().y
                + getOvYPosition() + getOvMaxHeight() - 1 ) );
    }

    final boolean inOvRectangle( final MouseEvent e ) {
        return ( ( e.getX() >= getOvRectangle().getX() - 1 )
                && ( e.getX() <= getOvRectangle().getX() + getOvRectangle().getWidth() + 1 )
                && ( e.getY() >= getOvRectangle().getY() - 1 ) && ( e.getY() <= getOvRectangle().getY()
                + getOvRectangle().getHeight() + 1 ) );
    }

    final private boolean inOvVirtualRectangle( final int x, final int y ) {
        return ( ( x >= getOvVirtualRectangle().x - 1 )
                && ( x <= getOvVirtualRectangle().x + getOvVirtualRectangle().width + 1 )
                && ( y >= getOvVirtualRectangle().y - 1 ) && ( y <= getOvVirtualRectangle().y
                + getOvVirtualRectangle().height + 1 ) );
    }

    final private boolean inOvVirtualRectangle( final MouseEvent e ) {
        return ( inOvVirtualRectangle( e.getX(), e.getY() ) );
    }

    final boolean isApplet() {
        return getMainPanel() instanceof MainPanelApplets;
    }

    final private boolean isCanBlast( final PhylogenyNode node ) {
        return ( node.getNodeData().isHasSequence() && ( ( ( node.getNodeData().getSequence().getAccession() != null ) && !ForesterUtil
                .isEmpty( node.getNodeData().getSequence().getAccession().getValue() ) )
                || !ForesterUtil.isEmpty( node.getNodeData().getSequence().getName() ) || !ForesterUtil.isEmpty( node
                .getNodeData().getSequence().getSymbol() ) ) );
    }

    final boolean isCanCollapse() {
        return ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
    }

    final boolean isCanColorSubtree() {
        return ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
    }

    final boolean isCanCopy() {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && getOptions().isEditable() );
    }

    final boolean isCanCut( final PhylogenyNode node ) {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && getOptions().isEditable() && !node
                .isRoot() );
    }

    final boolean isCanDelete() {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && getOptions().isEditable() );
    }

    final private boolean isCanOpenSeqWeb( final PhylogenyNode node ) {
        if ( node.getNodeData().isHasSequence()
                && ( node.getNodeData().getSequence().getAccession() != null )
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() )
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() )
                && getConfiguration().isHasWebLink( node.getNodeData().getSequence().getAccession().getSource()
                        .toLowerCase() ) ) {
            return true;
        }
        return false;
    }

    final private boolean isCanOpenTaxWeb( final PhylogenyNode node ) {
        if ( node.getNodeData().isHasTaxonomy()
                && ( ( ( node.getNodeData().getTaxonomy().getIdentifier() != null )
                        && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getIdentifier().getProvider() )
                        && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getIdentifier().getValue() ) && getConfiguration()
                        .isHasWebLink( node.getNodeData().getTaxonomy().getIdentifier().getProvider().toLowerCase() ) )
                        || ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) )
                        || ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getTaxonomyCode() ) ) || ( !ForesterUtil
                        .isEmpty( node.getNodeData().getTaxonomy().getCommonName() ) ) ) ) {
            return true;
        }
        else {
            return false;
        }
    }

    final boolean isCanPaste() {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && getOptions().isEditable()
                && ( getCutOrCopiedTree() != null ) && !getCutOrCopiedTree().isEmpty() );
    }

    final boolean isCanReroot() {
        return ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
    }

    final boolean isCanSubtree( final PhylogenyNode node ) {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && !node.isExternal() && ( !node
                .isRoot() || ( _subtree_index > 0 ) ) );
    }

    final boolean isEdited() {
        return _edited;
    }

    final private boolean isInFoundNodes( final PhylogenyNode node ) {
        return ( ( getFoundNodes() != null ) && getFoundNodes().contains( node ) );
    }

    final private boolean isInOv() {
        return _in_ov;
    }

    final boolean isInOvRect() {
        return _in_ov_rect;
    }

    final private boolean isNodeDataInvisible( final PhylogenyNode node ) {
        return ( ( node.getYcoord() < getVisibleRect().getMinY() - 40 )
                || ( node.getYcoord() > getVisibleRect().getMaxY() + 40 ) || ( ( node.getParent() != null ) && ( node
                .getParent().getXcoord() > getVisibleRect().getMaxX() ) ) );
    }

    final private boolean isNodeDataInvisibleUnrootedCirc( final PhylogenyNode node ) {
        return ( ( node.getYcoord() < getVisibleRect().getMinY() - 20 )
                || ( node.getYcoord() > getVisibleRect().getMaxY() + 20 )
                || ( node.getXcoord() < getVisibleRect().getMinX() - 20 ) || ( node.getXcoord() > getVisibleRect()
                .getMaxX() + 20 ) );
    }

    final private boolean isNonLinedUpCladogram() {
        return getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP;
    }

    final boolean isOvOn() {
        return _ov_on;
    }

    final boolean isPhyHasBranchLengths() {
        return _phy_has_branch_lengths;
    }

    final private boolean isUniformBranchLengthsForCladogram() {
        return getOptions().getCladogramType() == CLADOGRAM_TYPE.TOTAL_NODE_SUM_DEP;
    }

    final private void keyPressedCalls( final KeyEvent e ) {
        if ( isOvOn() && ( getMousePosition() != null ) && ( getMousePosition().getLocation() != null ) ) {
            if ( inOvVirtualRectangle( getMousePosition().x, getMousePosition().y ) ) {
                if ( !isInOvRect() ) {
                    setInOvRect( true );
                }
            }
            else if ( isInOvRect() ) {
                setInOvRect( false );
            }
        }
        if ( e.getModifiersEx() == InputEvent.CTRL_DOWN_MASK ) {
            if ( ( e.getKeyCode() == KeyEvent.VK_DELETE ) || ( e.getKeyCode() == KeyEvent.VK_HOME )
                    || ( e.getKeyCode() == KeyEvent.VK_F ) ) {
                getMainPanel().getTreeFontSet().mediumFonts();
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
            else if ( ( e.getKeyCode() == KeyEvent.VK_SUBTRACT ) || ( e.getKeyCode() == KeyEvent.VK_MINUS ) ) {
                getMainPanel().getTreeFontSet().decreaseFontSize();
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
            else if ( plusPressed( e.getKeyCode() ) ) {
                getMainPanel().getTreeFontSet().increaseFontSize();
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
        }
        else {
            if ( ( e.getKeyCode() == KeyEvent.VK_DELETE ) || ( e.getKeyCode() == KeyEvent.VK_HOME )
                    || ( e.getKeyCode() == KeyEvent.VK_F ) ) {
                getControlPanel().showWhole();
            }
            else if ( ( e.getKeyCode() == KeyEvent.VK_UP ) || ( e.getKeyCode() == KeyEvent.VK_DOWN )
                    || ( e.getKeyCode() == KeyEvent.VK_LEFT ) || ( e.getKeyCode() == KeyEvent.VK_RIGHT ) ) {
                if ( e.getModifiersEx() == InputEvent.SHIFT_DOWN_MASK ) {
                    if ( e.getKeyCode() == KeyEvent.VK_UP ) {
                        getMainPanel().getControlPanel().zoomInY( Constants.WHEEL_ZOOM_IN_FACTOR );
                        getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_DOWN ) {
                        getMainPanel().getControlPanel().zoomOutY( Constants.WHEEL_ZOOM_OUT_FACTOR );
                        getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_LEFT ) {
                        getMainPanel().getControlPanel().zoomOutX( Constants.WHEEL_ZOOM_OUT_FACTOR,
                                                                   Constants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR );
                        getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_RIGHT ) {
                        getMainPanel().getControlPanel().zoomInX( Constants.WHEEL_ZOOM_IN_FACTOR,
                                                                  Constants.WHEEL_ZOOM_IN_FACTOR );
                        getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
                else {
                    final int d = 80;
                    int dx = 0;
                    int dy = -d;
                    if ( e.getKeyCode() == KeyEvent.VK_DOWN ) {
                        dy = d;
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_LEFT ) {
                        dx = -d;
                        dy = 0;
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_RIGHT ) {
                        dx = d;
                        dy = 0;
                    }
                    final Point scroll_position = getMainPanel().getCurrentScrollPane().getViewport().getViewPosition();
                    scroll_position.x = scroll_position.x + dx;
                    scroll_position.y = scroll_position.y + dy;
                    if ( scroll_position.x <= 0 ) {
                        scroll_position.x = 0;
                    }
                    else {
                        final int max_x = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getMaximum()
                                - getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getVisibleAmount();
                        if ( scroll_position.x >= max_x ) {
                            scroll_position.x = max_x;
                        }
                    }
                    if ( scroll_position.y <= 0 ) {
                        scroll_position.y = 0;
                    }
                    else {
                        final int max_y = getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getMaximum()
                                - getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getVisibleAmount();
                        if ( scroll_position.y >= max_y ) {
                            scroll_position.y = max_y;
                        }
                    }
                    repaint();
                    getMainPanel().getCurrentScrollPane().getViewport().setViewPosition( scroll_position );
                }
            }
            else if ( ( e.getKeyCode() == KeyEvent.VK_SUBTRACT ) || ( e.getKeyCode() == KeyEvent.VK_MINUS ) ) {
                getMainPanel().getControlPanel().zoomOutY( Constants.WHEEL_ZOOM_OUT_FACTOR );
                getMainPanel().getControlPanel().zoomOutX( Constants.WHEEL_ZOOM_OUT_FACTOR,
                                                           Constants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR );
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
            }
            else if ( plusPressed( e.getKeyCode() ) ) {
                getMainPanel().getControlPanel().zoomInX( Constants.WHEEL_ZOOM_IN_FACTOR,
                                                          Constants.WHEEL_ZOOM_IN_FACTOR );
                getMainPanel().getControlPanel().zoomInY( Constants.WHEEL_ZOOM_IN_FACTOR );
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
            }
            else if ( e.getKeyCode() == KeyEvent.VK_S ) {
                if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                        || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                    setStartingAngle( ( getStartingAngle() % TWO_PI ) + ANGLE_ROTATION_UNIT );
                    getControlPanel().displayedPhylogenyMightHaveChanged( false );
                }
            }
            else if ( e.getKeyCode() == KeyEvent.VK_A ) {
                if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                        || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                    setStartingAngle( ( getStartingAngle() % TWO_PI ) - ANGLE_ROTATION_UNIT );
                    if ( getStartingAngle() < 0 ) {
                        setStartingAngle( TWO_PI + getStartingAngle() );
                    }
                    getControlPanel().displayedPhylogenyMightHaveChanged( false );
                }
            }
            else if ( e.getKeyCode() == KeyEvent.VK_D ) {
                boolean selected = false;
                if ( getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.HORIZONTAL ) {
                    getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
                    selected = true;
                }
                else {
                    getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.HORIZONTAL );
                }
                if ( getMainPanel().getMainFrame() == null ) {
                    // Must be "E" applet version.
                    final ArchaeopteryxE ae = ( ArchaeopteryxE ) ( ( MainPanelApplets ) getMainPanel() ).getApplet();
                    if ( ae.getlabelDirectionCbmi() != null ) {
                        ae.getlabelDirectionCbmi().setSelected( selected );
                    }
                }
                else {
                    getMainPanel().getMainFrame().getlabelDirectionCbmi().setSelected( selected );
                }
                repaint();
            }
            else if ( e.getKeyCode() == KeyEvent.VK_X ) {
                switchDisplaygetPhylogenyGraphicsType();
                repaint();
            }
            else if ( e.getKeyCode() == KeyEvent.VK_C ) {
                cycleColors();
                repaint();
            }
            else if ( getOptions().isShowOverview() && isOvOn() && ( e.getKeyCode() == KeyEvent.VK_O ) ) {
                MainFrame.cycleOverview( getOptions(), this );
                repaint();
            }
            else if ( getOptions().isShowOverview() && isOvOn() && ( e.getKeyCode() == KeyEvent.VK_I ) ) {
                increaseOvSize();
            }
            else if ( getOptions().isShowOverview() && isOvOn() && ( e.getKeyCode() == KeyEvent.VK_U ) ) {
                decreaseOvSize();
            }
            e.consume();
        }
    }

    final private void makePopupMenus( final PhylogenyNode node ) {
        _node_popup_menu = new JPopupMenu();
        final List<String> clickto_names = _main_panel.getControlPanel().getSingleClickToNames();
        _node_popup_menu_items = new JMenuItem[ clickto_names.size() ];
        for( int i = 0; i < clickto_names.size(); i++ ) {
            final String title = clickto_names.get( i );
            _node_popup_menu_items[ i ] = new JMenuItem( title );
            if ( title.equals( Configuration.clickto_options[ Configuration.open_seq_web ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanOpenSeqWeb( node ) );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.open_tax_web ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanOpenTaxWeb( node ) );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.blast ][ 0 ] ) ) {
                if ( Constants.__RELEASE || Constants.__SNAPSHOT_RELEASE ) {
                    continue;
                }
                _node_popup_menu_items[ i ].setEnabled( isCanBlast( node ) );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.delete_subtree_or_node ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
                _node_popup_menu_items[ i ].setEnabled( isCanDelete() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.cut_subtree ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
                _node_popup_menu_items[ i ].setEnabled( isCanCut( node ) );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.copy_subtree ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
                _node_popup_menu_items[ i ].setEnabled( isCanCopy() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.paste_subtree ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
                _node_popup_menu_items[ i ].setEnabled( isCanPaste() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.edit_node_data ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.add_new_node ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.reroot ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanReroot() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.collapse_uncollapse ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanCollapse() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.color_subtree ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanColorSubtree() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.subtree ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanSubtree( node ) );
            }
            _node_popup_menu_items[ i ].addActionListener( this );
            _node_popup_menu.add( _node_popup_menu_items[ i ] );
        }
    }

    final void midpointRoot() {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        if ( !_phylogeny.isRerootable() ) {
            JOptionPane.showMessageDialog( this,
                                           "This is not rerootable",
                                           "Not rerootable",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        setWaitCursor();
        PhylogenyMethods.midpointRoot( _phylogeny );
        resetNodeIdToDistToLeafMap();
        setArrowCursor();
        repaint();
    }

    final void mouseClicked( final MouseEvent e ) {
        if ( getOptions().isShowOverview() && isOvOn() && isInOv() ) {
            final double w_ratio = getVisibleRect().width / getOvRectangle().getWidth();
            final double h_ratio = getVisibleRect().height / getOvRectangle().getHeight();
            double x = ( e.getX() - getVisibleRect().x - getOvXPosition() - getOvRectangle().getWidth() / 2.0 )
                    * w_ratio;
            double y = ( e.getY() - getVisibleRect().y - getOvYPosition() - getOvRectangle().getHeight() / 2.0 )
                    * h_ratio;
            if ( x < 0 ) {
                x = 0;
            }
            if ( y < 0 ) {
                y = 0;
            }
            final double max_x = getWidth() - getVisibleRect().width;
            final double max_y = getHeight() - getVisibleRect().height;
            if ( x > max_x ) {
                x = max_x;
            }
            if ( y > max_y ) {
                y = max_y;
            }
            getMainPanel().getCurrentScrollPane().getViewport()
                    .setViewPosition( new Point( ForesterUtil.roundToInt( x ), ForesterUtil.roundToInt( y ) ) );
            setInOvRect( true );
            repaint();
        }
        else {
            final PhylogenyNode node = findNode( e.getX(), e.getY() );
            if ( node != null ) {
                if ( !node.isRoot() && node.getParent().isCollapse() ) {
                    return;
                }
                _highlight_node = node;
                // Check if shift key is down
                if ( ( e.getModifiers() & InputEvent.SHIFT_MASK ) != 0 ) {
                    // Yes, so add to _found_nodes
                    if ( getFoundNodes() == null ) {
                        setFoundNodes( new HashSet<PhylogenyNode>() );
                    }
                    getFoundNodes().add( node );
                    // Check if control key is down
                }
                else if ( ( e.getModifiers() & InputEvent.CTRL_MASK ) != 0 ) {
                    // Yes, so pop-up menu
                    displayNodePopupMenu( node, e.getX(), e.getY() );
                    // Handle unadorned click
                }
                else {
                    // Check for right mouse button
                    if ( e.getModifiers() == 4 ) {
                        displayNodePopupMenu( node, e.getX(), e.getY() );
                    }
                    else {
                        // if not in _found_nodes, clear _found_nodes
                        handleClickToAction( _control_panel.getActionWhenNodeClicked(), node );
                    }
                }
            }
            else {
                // no node was clicked
                _highlight_node = null;
            }
        }
        repaint();
    }

    final void mouseDragInBrowserPanel( final MouseEvent e ) {
        setCursor( MOVE_CURSOR );
        final Point scroll_position = getMainPanel().getCurrentScrollPane().getViewport().getViewPosition();
        scroll_position.x -= ( e.getX() - getLastDragPointX() );
        scroll_position.y -= ( e.getY() - getLastDragPointY() );
        if ( scroll_position.x < 0 ) {
            scroll_position.x = 0;
        }
        else {
            final int max_x = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getVisibleAmount();
            if ( scroll_position.x > max_x ) {
                scroll_position.x = max_x;
            }
        }
        if ( scroll_position.y < 0 ) {
            scroll_position.y = 0;
        }
        else {
            final int max_y = getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getVisibleAmount();
            if ( scroll_position.y > max_y ) {
                scroll_position.y = max_y;
            }
        }
        if ( isOvOn() || getOptions().isShowScale() ) {
            repaint();
        }
        getMainPanel().getCurrentScrollPane().getViewport().setViewPosition( scroll_position );
    }

    final void mouseDragInOvRectangle( final MouseEvent e ) {
        setCursor( HAND_CURSOR );
        final double w_ratio = getVisibleRect().width / getOvRectangle().getWidth();
        final double h_ratio = getVisibleRect().height / getOvRectangle().getHeight();
        final Point scroll_position = getMainPanel().getCurrentScrollPane().getViewport().getViewPosition();
        double dx = ( w_ratio * e.getX() - w_ratio * getLastDragPointX() );
        double dy = ( h_ratio * e.getY() - h_ratio * getLastDragPointY() );
        scroll_position.x = ForesterUtil.roundToInt( scroll_position.x + dx );
        scroll_position.y = ForesterUtil.roundToInt( scroll_position.y + dy );
        if ( scroll_position.x <= 0 ) {
            scroll_position.x = 0;
            dx = 0;
        }
        else {
            final int max_x = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getVisibleAmount();
            if ( scroll_position.x >= max_x ) {
                dx = 0;
                scroll_position.x = max_x;
            }
        }
        if ( scroll_position.y <= 0 ) {
            dy = 0;
            scroll_position.y = 0;
        }
        else {
            final int max_y = getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getVisibleAmount();
            if ( scroll_position.y >= max_y ) {
                dy = 0;
                scroll_position.y = max_y;
            }
        }
        repaint();
        getMainPanel().getCurrentScrollPane().getViewport().setViewPosition( scroll_position );
        setLastMouseDragPointX( ( float ) ( e.getX() + dx ) );
        setLastMouseDragPointY( ( float ) ( e.getY() + dy ) );
    }

    final void mouseMoved( final MouseEvent e ) {
        requestFocusInWindow();
        if ( getControlPanel().isNodeDescPopup() ) {
            if ( _node_desc_popup != null ) {
                _node_desc_popup.hide();
                _node_desc_popup = null;
            }
        }
        if ( getOptions().isShowOverview() && isOvOn() ) {
            if ( inOvVirtualRectangle( e ) ) {
                if ( !isInOvRect() ) {
                    setInOvRect( true );
                    repaint();
                }
            }
            else {
                if ( isInOvRect() ) {
                    setInOvRect( false );
                    repaint();
                }
            }
        }
        if ( inOv( e ) && getOptions().isShowOverview() && isOvOn() ) {
            if ( !isInOv() ) {
                setInOv( true );
            }
        }
        else {
            if ( isInOv() ) {
                setInOv( false );
            }
            final PhylogenyNode node = findNode( e.getX(), e.getY() );
            if ( ( node != null ) && ( node.isRoot() || !node.getParent().isCollapse() ) ) {
                // cursor is over a tree node
                if ( ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.CUT_SUBTREE )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.COPY_SUBTREE )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.PASTE_SUBTREE )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.DELETE_NODE_OR_SUBTREE )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.REROOT )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.ADD_NEW_NODE ) ) {
                    setCursor( CUT_CURSOR );
                }
                else {
                    setCursor( HAND_CURSOR );
                    if ( getControlPanel().isNodeDescPopup() ) {
                        showNodeDataPopup( e, node );
                    }
                }
            }
            else {
                setCursor( ARROW_CURSOR );
            }
        }
    }

    final void mouseReleasedInBrowserPanel( final MouseEvent e ) {
        setCursor( ARROW_CURSOR );
    }

    final public void mouseWheelMoved( final MouseWheelEvent e ) {
        final int notches = e.getWheelRotation();
        if ( inOvVirtualRectangle( e ) ) {
            if ( !isInOvRect() ) {
                setInOvRect( true );
                repaint();
            }
        }
        else {
            if ( isInOvRect() ) {
                setInOvRect( false );
                repaint();
            }
        }
        if ( e.isControlDown() ) {
            if ( notches < 0 ) {
                getTreeFontSet().increaseFontSize();
                getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
            else {
                getTreeFontSet().decreaseFontSize();
                getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
        }
        else if ( e.isShiftDown() ) {
            if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                    || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                if ( notches < 0 ) {
                    for( int i = 0; i < ( -notches ); ++i ) {
                        setStartingAngle( ( getStartingAngle() % TWO_PI ) + ANGLE_ROTATION_UNIT );
                        getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
                else {
                    for( int i = 0; i < notches; ++i ) {
                        setStartingAngle( ( getStartingAngle() % TWO_PI ) - ANGLE_ROTATION_UNIT );
                        if ( getStartingAngle() < 0 ) {
                            setStartingAngle( TWO_PI + getStartingAngle() );
                        }
                        getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
            }
            else {
                if ( notches < 0 ) {
                    for( int i = 0; i < ( -notches ); ++i ) {
                        getControlPanel().zoomInY( Constants.WHEEL_ZOOM_IN_FACTOR );
                        getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
                else {
                    for( int i = 0; i < notches; ++i ) {
                        getControlPanel().zoomOutY( Constants.WHEEL_ZOOM_OUT_FACTOR );
                        getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
            }
        }
        else {
            if ( notches < 0 ) {
                for( int i = 0; i < ( -notches ); ++i ) {
                    getControlPanel().zoomInX( Constants.WHEEL_ZOOM_IN_FACTOR,
                                               Constants.WHEEL_ZOOM_IN_X_CORRECTION_FACTOR );
                    getControlPanel().zoomInY( Constants.WHEEL_ZOOM_IN_FACTOR );
                    getControlPanel().displayedPhylogenyMightHaveChanged( false );
                }
            }
            else {
                for( int i = 0; i < notches; ++i ) {
                    getControlPanel().zoomOutY( Constants.WHEEL_ZOOM_OUT_FACTOR );
                    getControlPanel().zoomOutX( Constants.WHEEL_ZOOM_OUT_FACTOR,
                                                Constants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR );
                    getControlPanel().displayedPhylogenyMightHaveChanged( false );
                }
            }
        }
        requestFocus();
        requestFocusInWindow();
        requestFocus();
    }

    final void multiplyUrtFactor( final float f ) {
        _urt_factor *= f;
    }

    final JApplet obtainApplet() {
        return ( ( MainPanelApplets ) getMainPanel() ).getApplet();
    }

    final private void openSeqWeb( final PhylogenyNode node ) {
        if ( !isCanOpenSeqWeb( node ) ) {
            cannotOpenBrowserWarningMessage( "sequence" );
            return;
        }
        String uri_str = null;
        final Sequence seq = node.getNodeData().getSequence();
        final String source = seq.getAccession().getSource().toLowerCase();
        final WebLink weblink = getConfiguration().getWebLink( source );
        try {
            uri_str = weblink.getUrl() + URLEncoder.encode( seq.getAccession().getValue(), ForesterConstants.UTF8 );
        }
        catch ( final UnsupportedEncodingException e ) {
            Util.showErrorMessage( this, e.toString() );
            e.printStackTrace();
        }
        if ( !ForesterUtil.isEmpty( uri_str ) ) {
            try {
                JApplet applet = null;
                if ( isApplet() ) {
                    applet = obtainApplet();
                }
                Util.launchWebBrowser( new URI( uri_str ), isApplet(), applet, "_aptx_seq" );
            }
            catch ( final IOException e ) {
                Util.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
            catch ( final URISyntaxException e ) {
                Util.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        else {
            cannotOpenBrowserWarningMessage( "sequence" );
        }
    }

    final private void openTaxWeb( final PhylogenyNode node ) {
        if ( !isCanOpenTaxWeb( node ) ) {
            cannotOpenBrowserWarningMessage( "taxonomic" );
            return;
        }
        String uri_str = null;
        final Taxonomy tax = node.getNodeData().getTaxonomy();
        if ( ( tax.getIdentifier() != null ) && !ForesterUtil.isEmpty( tax.getIdentifier().getProvider() ) ) {
            final String type = tax.getIdentifier().getProvider().toLowerCase();
            if ( getConfiguration().isHasWebLink( type ) ) {
                final WebLink weblink = getConfiguration().getWebLink( type );
                try {
                    uri_str = weblink.getUrl()
                            + URLEncoder.encode( tax.getIdentifier().getValue(), ForesterConstants.UTF8 );
                }
                catch ( final UnsupportedEncodingException e ) {
                    Util.showErrorMessage( this, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        else if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            try {
                uri_str = "http://www.eol.org/search?q="
                        + URLEncoder.encode( tax.getScientificName(), ForesterConstants.UTF8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                Util.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        else if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
            try {
                uri_str = "http://www.uniprot.org/taxonomy/?query="
                        + URLEncoder.encode( tax.getTaxonomyCode(), ForesterConstants.UTF8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                Util.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        else if ( !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
            try {
                uri_str = "http://www.eol.org/search?q="
                        + URLEncoder.encode( tax.getCommonName(), ForesterConstants.UTF8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                Util.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        if ( !ForesterUtil.isEmpty( uri_str ) ) {
            try {
                JApplet applet = null;
                if ( isApplet() ) {
                    applet = obtainApplet();
                }
                Util.launchWebBrowser( new URI( uri_str ), isApplet(), applet, "_aptx_tax" );
            }
            catch ( final IOException e ) {
                Util.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
            catch ( final URISyntaxException e ) {
                Util.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        else {
            cannotOpenBrowserWarningMessage( "taxonomic" );
        }
    }

    final void paintBranchCircular( final PhylogenyNode p,
                                    final PhylogenyNode c,
                                    final Graphics2D g,
                                    final boolean radial_labels,
                                    final boolean to_pdf,
                                    final boolean to_graphics_file ) {
        final double angle = _urt_nodeid_angle_map.get( c.getNodeId() );
        final double root_x = _root.getXcoord();
        final double root_y = _root.getYcoord();
        final double dx = root_x - p.getXcoord();
        final double dy = root_y - p.getYcoord();
        final double parent_radius = Math.sqrt( dx * dx + dy * dy );
        final double arc = ( _urt_nodeid_angle_map.get( p.getNodeId() ) ) - angle;
        assignGraphicsForBranchWithColorForParentBranch( c, false, g, to_pdf, to_graphics_file );
        if ( ( c.isFirstChildNode() || c.isLastChildNode() )
                && ( ( Math.abs( parent_radius * arc ) > 1.5 ) || to_pdf || to_graphics_file ) ) {
            final double r2 = 2.0 * parent_radius;
            drawArc( root_x - parent_radius, root_y - parent_radius, r2, r2, ( -angle - arc ), arc, g );
        }
        drawLine( c.getXcoord(), c.getYcoord(), root_x + ( Math.cos( angle ) * parent_radius ), root_y
                + ( Math.sin( angle ) * parent_radius ), g );
        paintNodeBox( c.getXcoord(), c.getYcoord(), c, g, to_pdf, to_graphics_file, isInFoundNodes( c ) );
        if ( c.isExternal() ) {
            final boolean is_in_found_nodes = isInFoundNodes( c );
            if ( ( _dynamic_hiding_factor > 1 ) && !is_in_found_nodes
                    && ( _urt_nodeid_index_map.get( c.getNodeId() ) % _dynamic_hiding_factor != 1 ) ) {
                return;
            }
            paintNodeDataUnrootedCirc( g, c, to_pdf, to_graphics_file, radial_labels, 0, is_in_found_nodes );
        }
    }

    final void paintBranchCircularLite( final PhylogenyNode p, final PhylogenyNode c, final Graphics2D g ) {
        final double angle = _urt_nodeid_angle_map.get( c.getNodeId() );
        final double root_x = _root.getXSecondary();
        final double root_y = _root.getYSecondary();
        final double dx = root_x - p.getXSecondary();
        final double dy = root_y - p.getYSecondary();
        final double arc = ( _urt_nodeid_angle_map.get( p.getNodeId() ) ) - angle;
        final double parent_radius = Math.sqrt( dx * dx + dy * dy );
        g.setColor( getTreeColorSet().getOvColor() );
        if ( ( c.isFirstChildNode() || c.isLastChildNode() ) && ( Math.abs( arc ) > 0.02 ) ) {
            final double r2 = 2.0 * parent_radius;
            drawArc( root_x - parent_radius, root_y - parent_radius, r2, r2, ( -angle - arc ), arc, g );
        }
        drawLine( c.getXSecondary(), c.getYSecondary(), root_x + ( Math.cos( angle ) * parent_radius ), root_y
                + ( Math.sin( angle ) * parent_radius ), g );
        if ( isInFoundNodes( c ) ) {
            g.setColor( getTreeColorSet().getFoundColor() );
            drawRectFilled( c.getXSecondary() - 1, c.getYSecondary() - 1, 3, 3, g );
        }
    }

    final private void paintBranchLength( final Graphics2D g,
                                          final PhylogenyNode node,
                                          final boolean to_pdf,
                                          final boolean to_graphics_file ) {
        g.setFont( getTreeFontSet().getSmallFont() );
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else {
            g.setColor( getTreeColorSet().getBranchLengthColor() );
        }
        if ( !node.isRoot() ) {
            if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                TreePanel.drawString( FORMATTER_BRANCH_LENGTH.format( node.getDistanceToParent() ), node.getParent()
                        .getXcoord()
                        + EURO_D, node.getYcoord() - getTreeFontSet()._small_max_descent, g );
            }
            else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                TreePanel.drawString( FORMATTER_BRANCH_LENGTH.format( node.getDistanceToParent() ), node.getParent()
                        .getXcoord()
                        + ROUNDED_D, node.getYcoord() - getTreeFontSet()._small_max_descent, g );
            }
            else {
                TreePanel.drawString( FORMATTER_BRANCH_LENGTH.format( node.getDistanceToParent() ), node.getParent()
                        .getXcoord() + 3, node.getYcoord() - getTreeFontSet()._small_max_descent, g );
            }
        }
        else {
            TreePanel.drawString( FORMATTER_BRANCH_LENGTH.format( node.getDistanceToParent() ), 3, node.getYcoord()
                    - getTreeFontSet()._small_max_descent, g );
        }
    }

    final private void paintBranchLite( final Graphics2D g,
                                        final float x1,
                                        final float x2,
                                        final float y1,
                                        final float y2,
                                        final PhylogenyNode node ) {
        g.setColor( getTreeColorSet().getOvColor() );
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR ) {
            drawLine( x1, y1, x2, y2, g );
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CONVEX ) {
            _quad_curve.setCurve( x1, y1, x1, y2, x2, y2 );
            ( g ).draw( _quad_curve );
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CURVED ) {
            final float dx = x2 - x1;
            final float dy = y2 - y1;
            _cubic_curve.setCurve( x1, y1, x1 + ( dx * 0.4f ), y1 + ( dy * 0.2f ), x1 + ( dx * 0.6f ), y1
                    + ( dy * 0.8f ), x2, y2 );
            ( g ).draw( _cubic_curve );
        }
        else {
            final float x2a = x2;
            final float x1a = x1;
            // draw the vertical line
            if ( node.isFirstChildNode() || node.isLastChildNode() ) {
                drawLine( x1, y1, x1, y2, g );
            }
            // draw the horizontal line
            drawLine( x1a, y2, x2a, y2, g );
        }
    }

    /**
     * Paint a branch which consists of a vertical and a horizontal bar
     * @param is_ind_found_nodes 
     */
    final private void paintBranchRectangular( final Graphics2D g,
                                               final float x1,
                                               final float x2,
                                               float y1,
                                               final float y2,
                                               final PhylogenyNode node,
                                               final boolean to_pdf,
                                               final boolean to_graphics_file ) {
        assignGraphicsForBranchWithColorForParentBranch( node, false, g, to_pdf, to_graphics_file );
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR ) {
            drawLine( x1, y1, x2, y2, g );
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CONVEX ) {
            _quad_curve.setCurve( x1, y1, x1, y2, x2, y2 );
            g.draw( _quad_curve );
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CURVED ) {
            final float dx = x2 - x1;
            final float dy = y2 - y1;
            _cubic_curve.setCurve( x1, y1, x1 + ( dx * 0.4f ), y1 + ( dy * 0.2f ), x1 + ( dx * 0.6f ), y1
                    + ( dy * 0.8f ), x2, y2 );
            g.draw( _cubic_curve );
        }
        else {
            float x2a = x2;
            float x1a = x1;
            // draw the vertical line
            boolean draw_horizontal = true;
            float y2_r = 0;
            if ( node.isFirstChildNode() || node.isLastChildNode()
                    || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE )
                    || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) ) {
                boolean draw_vertical = true;
                final PhylogenyNode parent = node.getParent();
                if ( ( ( getOptions().isShowNodeBoxes() && !to_pdf && !to_graphics_file ) || ( ( getControlPanel()
                        .isEvents() )
                        && ( parent != null ) && parent.isHasAssignedEvent() ) )
                        && ( _phylogeny.isRooted() || !( ( parent != null ) && parent.isRoot() ) )
                        && !( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() && !parent
                                .isDuplication() ) ) {
                    if ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE )
                            && ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) ) {
                        if ( Math.abs( y2 - y1 ) <= TreePanel.HALF_BOX_SIZE ) {
                            draw_vertical = false;
                        }
                        else {
                            if ( y1 < y2 ) {
                                y1 += TreePanel.HALF_BOX_SIZE;
                            }
                            else {
                                if ( !to_pdf ) {
                                    y1 -= TreePanel.HALF_BOX_SIZE + 1;
                                }
                                else {
                                    y1 -= TreePanel.HALF_BOX_SIZE;
                                }
                            }
                        }
                    }
                    if ( ( x2 - x1 ) <= TreePanel.HALF_BOX_SIZE ) {
                        draw_horizontal = false;
                    }
                    else if ( !draw_vertical ) {
                        x1a += TreePanel.HALF_BOX_SIZE;
                    }
                    if ( ( ( x2 - x1a ) > TreePanel.HALF_BOX_SIZE )
                            && !( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() && !node
                                    .isDuplication() ) ) {
                        x2a -= TreePanel.HALF_BOX_SIZE;
                    }
                }
                if ( draw_vertical ) {
                    if ( !to_graphics_file
                            && !to_pdf
                            && ( ( ( y2 < getVisibleRect().getMinY() - 20 ) && ( y1 < getVisibleRect().getMinY() - 20 ) ) || ( ( y2 > getVisibleRect()
                                    .getMaxY() + 20 ) && ( y1 > getVisibleRect().getMaxY() + 20 ) ) ) ) {
                        // Do nothing.
                    }
                    else {
                        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                            float x2c = x1 + EURO_D;
                            if ( x2c > x2a ) {
                                x2c = x2a;
                            }
                            drawLine( x1, y1, x2c, y2, g );
                        }
                        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                            if ( y2 > y1 ) {
                                y2_r = y2 - ROUNDED_D;
                                if ( y2_r < y1 ) {
                                    y2_r = y1;
                                }
                                drawLine( x1, y1, x1, y2_r, g );
                            }
                            else {
                                y2_r = y2 + ROUNDED_D;
                                if ( y2_r > y1 ) {
                                    y2_r = y1;
                                }
                                drawLine( x1, y1, x1, y2_r, g );
                            }
                        }
                        else {
                            drawLine( x1, y1, x1, y2, g );
                        }
                    }
                }
            }
            // draw the horizontal line
            if ( !to_graphics_file && !to_pdf
                    && ( ( y2 < getVisibleRect().getMinY() - 20 ) || ( y2 > getVisibleRect().getMaxY() + 20 ) ) ) {
                return;
            }
            float x1_r = 0;
            if ( draw_horizontal ) {
                if ( !getControlPanel().isWidthBranches() || ( PhylogenyMethods.getBranchWidthValue( node ) == 1 ) ) {
                    if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                        x1_r = x1a + ROUNDED_D;
                        if ( x1_r < x2a ) {
                            drawLine( x1_r, y2, x2a, y2, g );
                        }
                    }
                    else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                        final float x1c = x1a + EURO_D;
                        if ( x1c < x2a ) {
                            drawLine( x1c, y2, x2a, y2, g );
                        }
                    }
                    else {
                        drawLine( x1a, y2, x2a, y2, g );
                    }
                }
                else {
                    final double w = PhylogenyMethods.getBranchWidthValue( node );
                    if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                        x1_r = x1a + ROUNDED_D;
                        if ( x1_r < x2a ) {
                            drawRectFilled( x1_r, y2 - ( w / 2 ), x2a - x1_r, w, g );
                        }
                    }
                    else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                        final float x1c = x1a + EURO_D;
                        if ( x1c < x2a ) {
                            drawRectFilled( x1c, y2 - ( w / 2 ), x2a - x1c, w, g );
                        }
                    }
                    else {
                        drawRectFilled( x1a, y2 - ( w / 2 ), x2a - x1a, w, g );
                    }
                }
            }
            if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) ) {
                if ( x1_r > x2a ) {
                    x1_r = x2a;
                }
                if ( y2 > y2_r ) {
                    final double diff = y2 - y2_r;
                    _arc.setArc( x1, y2_r - diff, 2 * ( x1_r - x1 ), 2 * diff, 180, 90, Arc2D.OPEN );
                }
                else {
                    _arc.setArc( x1, y2, 2 * ( x1_r - x1 ), 2 * ( y2_r - y2 ), 90, 90, Arc2D.OPEN );
                }
                g.draw( _arc );
            }
        }
        paintNodeBox( x2, y2, node, g, to_pdf, to_graphics_file, isInFoundNodes( node ) );
    }

    final void paintCircular( final Phylogeny phy,
                              final double starting_angle,
                              final int center_x,
                              final int center_y,
                              final int radius,
                              final Graphics2D g,
                              final boolean to_pdf,
                              final boolean to_graphics_file ) {
        _circ_num_ext_nodes = phy.getNumberOfExternalNodes();
        _root = phy.getRoot();
        _root.setXcoord( center_x );
        _root.setYcoord( center_y );
        paintNodeBox( _root.getXcoord(), _root.getYcoord(), _root, g, to_pdf, to_graphics_file, isInFoundNodes( _root ) );
        final boolean radial_labels = getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL;
        double current_angle = starting_angle;
        int i = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.setXcoord( ( float ) ( center_x + ( radius * Math.cos( current_angle ) ) ) );
            n.setYcoord( ( float ) ( center_y + ( radius * Math.sin( current_angle ) ) ) );
            _urt_nodeid_angle_map.put( n.getNodeId(), current_angle );
            _urt_nodeid_index_map.put( n.getNodeId(), i++ );
            current_angle += ( TWO_PI / _circ_num_ext_nodes );
        }
        paintCirculars( phy.getRoot(), phy, center_x, center_y, radius, radial_labels, g, to_pdf, to_graphics_file );
    }

    final void paintCircularLite( final Phylogeny phy,
                                  final double starting_angle,
                                  final int center_x,
                                  final int center_y,
                                  final int radius,
                                  final Graphics2D g ) {
        _circ_num_ext_nodes = phy.getNumberOfExternalNodes();
        _root = phy.getRoot();
        _root.setXSecondary( center_x );
        _root.setYSecondary( center_y );
        double current_angle = starting_angle;
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.setXSecondary( ( float ) ( center_x + radius * Math.cos( current_angle ) ) );
            n.setYSecondary( ( float ) ( center_y + radius * Math.sin( current_angle ) ) );
            _urt_nodeid_angle_map.put( n.getNodeId(), current_angle );
            current_angle += ( TWO_PI / _circ_num_ext_nodes );
        }
        paintCircularsLite( phy.getRoot(), phy, center_x, center_y, radius, g );
    }

    final private double paintCirculars( final PhylogenyNode n,
                                         final Phylogeny phy,
                                         final float center_x,
                                         final float center_y,
                                         final double radius,
                                         final boolean radial_labels,
                                         final Graphics2D g,
                                         final boolean to_pdf,
                                         final boolean to_graphics_file ) {
        if ( n.isExternal() ) {
            return _urt_nodeid_angle_map.get( n.getNodeId() );
        }
        else {
            final List<PhylogenyNode> descs = n.getDescendants();
            double sum = 0;
            for( final PhylogenyNode desc : descs ) {
                sum += paintCirculars( desc,
                                       phy,
                                       center_x,
                                       center_y,
                                       radius,
                                       radial_labels,
                                       g,
                                       to_pdf,
                                       to_graphics_file );
            }
            double r = 0;
            if ( !n.isRoot() ) {
                r = 1 - ( ( ( double ) _circ_max_depth - PhylogenyMethods.calculateDepth( n ) ) / _circ_max_depth );
            }
            final double theta = sum / descs.size();
            n.setXcoord( ( float ) ( center_x + r * radius * Math.cos( theta ) ) );
            n.setYcoord( ( float ) ( center_y + r * radius * Math.sin( theta ) ) );
            _urt_nodeid_angle_map.put( n.getNodeId(), theta );
            for( final PhylogenyNode desc : descs ) {
                paintBranchCircular( n, desc, g, radial_labels, to_pdf, to_graphics_file );
            }
            return theta;
        }
    }

    final private void paintCircularsLite( final PhylogenyNode n,
                                           final Phylogeny phy,
                                           final int center_x,
                                           final int center_y,
                                           final int radius,
                                           final Graphics2D g ) {
        if ( n.isExternal() ) {
            return;
        }
        else {
            final List<PhylogenyNode> descs = n.getDescendants();
            for( final PhylogenyNode desc : descs ) {
                paintCircularsLite( desc, phy, center_x, center_y, radius, g );
            }
            float r = 0;
            if ( !n.isRoot() ) {
                r = 1 - ( ( ( float ) _circ_max_depth - PhylogenyMethods.calculateDepth( n ) ) / _circ_max_depth );
            }
            final double theta = _urt_nodeid_angle_map.get( n.getNodeId() );
            n.setXSecondary( ( float ) ( center_x + radius * r * Math.cos( theta ) ) );
            n.setYSecondary( ( float ) ( center_y + radius * r * Math.sin( theta ) ) );
            for( final PhylogenyNode desc : descs ) {
                paintBranchCircularLite( n, desc, g );
            }
        }
    }

    final private void paintCollapsedNode( final Graphics2D g,
                                           final PhylogenyNode node,
                                           final boolean to_graphics_file,
                                           final boolean to_pdf,
                                           final boolean is_in_found_nodes ) {
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else if ( is_in_found_nodes ) {
            g.setColor( getTreeColorSet().getFoundColor() );
        }
        else if ( getControlPanel().isColorAccordingToTaxonomy() ) {
            g.setColor( getTaxonomyBasedColor( node ) );
        }
        else {
            g.setColor( getTreeColorSet().getCollapseFillColor() );
        }
        double d = node.getAllExternalDescendants().size();
        if ( d > 1000 ) {
            d = ( 3 * _y_distance ) / 3;
        }
        else {
            d = ( Math.log10( d ) * _y_distance ) / 2.5;
        }
        if ( d < BOX_SIZE ) {
            d = BOX_SIZE;
        }
        _polygon.reset();
        _polygon.addPoint( ForesterUtil.roundToInt( node.getXcoord() - TreePanel.BOX_SIZE ), ForesterUtil
                .roundToInt( node.getYcoord() ) );
        _polygon.addPoint( ForesterUtil.roundToInt( node.getXcoord() + TreePanel.BOX_SIZE ), ForesterUtil
                .roundToInt( node.getYcoord() - d ) );
        _polygon.addPoint( ForesterUtil.roundToInt( node.getXcoord() + TreePanel.BOX_SIZE ), ForesterUtil
                .roundToInt( node.getYcoord() + d ) );
        g.fillPolygon( _polygon );
        paintNodeData( g, node, to_graphics_file, to_pdf, is_in_found_nodes );
    }

    @Override
    final public void paintComponent( final Graphics g ) {
        super.paintComponent( g );
        final Graphics2D g2d = ( Graphics2D ) g;
        g2d.setRenderingHints( _rendering_hints );
        paintPhylogeny( g2d, false, false, 0, 0, 0, 0 );
    }

    final private void paintConfidenceValues( final Graphics2D g,
                                              final PhylogenyNode node,
                                              final boolean to_pdf,
                                              final boolean to_graphics_file ) {
        String conf_str = "";
        final List<Confidence> confidences = node.getBranchData().getConfidences();
        if ( confidences.size() == 1 ) {
            final double value = node.getBranchData().getConfidence( 0 ).getValue();
            if ( ( value == Confidence.CONFIDENCE_DEFAULT_VALUE ) || ( value < getOptions().getMinConfidenceValue() ) ) {
                return;
            }
            conf_str = FORMATTER_CONFIDENCE.format( value );
        }
        else if ( confidences.size() > 1 ) {
            boolean one_ok = false;
            boolean not_first = false;
            Collections.sort( confidences );
            final StringBuilder sb = new StringBuilder();
            for( final Confidence confidence : confidences ) {
                final double value = confidence.getValue();
                if ( value != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                    if ( value >= getOptions().getMinConfidenceValue() ) {
                        one_ok = true;
                    }
                    if ( not_first ) {
                        sb.append( "/" );
                    }
                    else {
                        not_first = true;
                    }
                    sb.append( FORMATTER_CONFIDENCE.format( ForesterUtil.round( value, getOptions()
                            .getNumberOfDigitsAfterCommaForConfidenceValues() ) ) );
                }
            }
            if ( one_ok ) {
                conf_str = sb.toString();
            }
        }
        if ( conf_str.length() > 0 ) {
            final double parent_x = node.getParent().getXcoord();
            double x = node.getXcoord();
            g.setFont( getTreeFontSet().getSmallFont() );
            if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                x += EURO_D;
            }
            else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                x += ROUNDED_D;
            }
            if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
                g.setColor( Color.BLACK );
            }
            else {
                g.setColor( getTreeColorSet().getConfidenceColor() );
            }
            TreePanel
                    .drawString( conf_str, parent_x
                            + ( ( x - parent_x - getTreeFontSet()._fm_small.stringWidth( conf_str ) ) / 2 ), ( node
                            .getYcoord() + getTreeFontSet()._small_max_ascent ) - 1, g );
        }
    }

    final private void paintFoundNode( final int x, final int y, final Graphics2D g ) {
        g.setColor( getTreeColorSet().getFoundColor() );
        g.fillRect( x - TreePanel.HALF_BOX_SIZE, y - TreePanel.HALF_BOX_SIZE, TreePanel.BOX_SIZE, TreePanel.BOX_SIZE );
    }

    final private void paintGainedAndLostCharacters( final Graphics2D g,
                                                     final PhylogenyNode node,
                                                     final String gained,
                                                     final String lost ) {
        if ( node.getParent() != null ) {
            final double parent_x = node.getParent().getXcoord();
            final double x = node.getXcoord();
            g.setFont( getTreeFontSet().getLargeFont() );
            g.setColor( getTreeColorSet().getGainedCharactersColor() );
            if ( Constants.SPECIAL_CUSTOM ) {
                g.setColor( Color.BLUE );
            }
            TreePanel
                    .drawString( gained, parent_x
                            + ( ( x - parent_x - getTreeFontSet()._fm_large.stringWidth( gained ) ) / 2 ), ( node
                            .getYcoord() - getTreeFontSet()._fm_large.getMaxDescent() ) - 1, g );
            g.setColor( getTreeColorSet().getLostCharactersColor() );
            TreePanel.drawString( lost,
                                  parent_x + ( ( x - parent_x - getTreeFontSet()._fm_large.stringWidth( lost ) ) / 2 ),
                                  ( node.getYcoord() + getTreeFontSet()._fm_large.getMaxAscent() ) + 1,
                                  g );
        }
    }

    /**
     * Draw a box at the indicated node.
     * 
     * @param x
     * @param y
     * @param node
     * @param g
     */
    final private void paintNodeBox( final double x,
                                     final double y,
                                     final PhylogenyNode node,
                                     final Graphics2D g,
                                     final boolean to_pdf,
                                     final boolean to_graphics_file,
                                     final boolean is_in_found_nodes ) {
        if ( node.isCollapse() ) {
            return;
        }
        // if this node should be highlighted, do so
        if ( ( _highlight_node == node ) && !to_pdf && !to_graphics_file ) {
            g.setColor( getTreeColorSet().getFoundColor() );
            drawOval( x - 8, y - 8, 16, 16, g );
            drawOval( x - 9, y - 8, 17, 17, g );
            drawOval( x - 9, y - 9, 18, 18, g );
        }
        if ( is_in_found_nodes ) {
            paintFoundNode( ForesterUtil.roundToInt( x ), ForesterUtil.roundToInt( y ), g );
        }
        else {
            if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
                g.setColor( Color.BLACK );
            }
            else if ( getControlPanel().isEvents() && Util.isHasAssignedEvent( node ) ) {
                final Event event = node.getNodeData().getEvent();
                if ( event.isDuplication() ) {
                    g.setColor( getTreeColorSet().getDuplicationBoxColor() );
                }
                else if ( event.isSpeciation() ) {
                    g.setColor( getTreeColorSet().getSpecBoxColor() );
                }
                else if ( event.isSpeciationOrDuplication() ) {
                    g.setColor( getTreeColorSet().getDuplicationOrSpeciationColor() );
                }
            }
            else {
                assignGraphicsForNodeBoxWithColorForParentBranch( node, g );
            }
            if ( ( getOptions().isShowNodeBoxes() && !to_pdf && !to_graphics_file )
                    || ( getControlPanel().isEvents() && node.isHasAssignedEvent() ) ) {
                if ( to_pdf || to_graphics_file ) {
                    if ( node.isDuplication() || !getOptions().isPrintBlackAndWhite() ) {
                        drawOvalFilled( x - HALF_BOX_SIZE, y - HALF_BOX_SIZE, BOX_SIZE, BOX_SIZE, g );
                    }
                }
                else {
                    drawRectFilled( x - HALF_BOX_SIZE, y - HALF_BOX_SIZE, BOX_SIZE, BOX_SIZE, g );
                }
            }
        }
    }

    final private void paintNodeData( final Graphics2D g,
                                      final PhylogenyNode node,
                                      final boolean to_graphics_file,
                                      final boolean to_pdf,
                                      final boolean is_in_found_nodes ) {
        if ( isNodeDataInvisible( node ) && !to_graphics_file && !to_pdf ) {
            return;
        }
        if ( getOptions().isShowBranchLengthValues()
                && ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR )
                        || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) )
                && ( !node.isRoot() ) && ( node.getDistanceToParent() != PhylogenyNode.DISTANCE_DEFAULT ) ) {
            paintBranchLength( g, node, to_pdf, to_graphics_file );
        }
        if ( !getControlPanel().isShowInternalData() && !node.isExternal() && !node.isCollapse() ) {
            return;
        }
        int x = 0;
        if ( node.getNodeData().isHasTaxonomy()
                && ( getControlPanel().isShowTaxonomyCode() || getControlPanel().isShowTaxonomyNames() ) ) {
            x = paintTaxonomy( g, node, is_in_found_nodes, to_pdf, to_graphics_file );
        }
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else if ( is_in_found_nodes ) {
            g.setColor( getTreeColorSet().getFoundColor() );
        }
        else if ( getControlPanel().isColorAccordingToTaxonomy() ) {
            g.setColor( getTaxonomyBasedColor( node ) );
        }
        else {
            g.setColor( getTreeColorSet().getSequenceColor() );
        }
        _sb.setLength( 0 );
        if ( node.isCollapse() && ( ( !node.isRoot() && !node.getParent().isCollapse() ) || node.isRoot() ) ) {
            _sb.append( " [" );
            _sb.append( node.getAllExternalDescendants().size() );
            _sb.append( "]" );
        }
        if ( getControlPanel().isShowNodeNames() && ( node.getNodeName().length() > 0 ) ) {
            if ( _sb.length() > 0 ) {
                _sb.append( " " );
            }
            _sb.append( node.getNodeName() );
        }
        if ( node.getNodeData().isHasSequence() ) {
            if ( getControlPanel().isShowGeneSymbols() && ( node.getNodeData().getSequence().getSymbol().length() > 0 ) ) {
                if ( _sb.length() > 0 ) {
                    _sb.append( " " );
                }
                _sb.append( node.getNodeData().getSequence().getSymbol() );
            }
            if ( getControlPanel().isShowGeneNames() && ( node.getNodeData().getSequence().getName().length() > 0 ) ) {
                if ( _sb.length() > 0 ) {
                    _sb.append( " " );
                }
                _sb.append( node.getNodeData().getSequence().getName() );
            }
            if ( getControlPanel().isShowSequenceAcc() && ( node.getNodeData().getSequence().getAccession() != null ) ) {
                if ( _sb.length() > 0 ) {
                    _sb.append( " " );
                }
                if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() ) ) {
                    _sb.append( node.getNodeData().getSequence().getAccession().getSource() );
                    _sb.append( ":" );
                }
                _sb.append( node.getNodeData().getSequence().getAccession().getValue() );
            }
        }
        g.setFont( getTreeFontSet().getLargeFont() );
        if ( is_in_found_nodes ) {
            g.setFont( getTreeFontSet().getLargeFont().deriveFont( Font.BOLD ) );
        }
        double down_shift_factor = 3.0;
        if ( !node.isExternal() && ( node.getNumberOfDescendants() == 1 ) ) {
            down_shift_factor = 1;
        }
        if ( _sb.length() > 0 ) {
            TreePanel.drawString( _sb.toString(), node.getXcoord() + x + 2 + TreePanel.HALF_BOX_SIZE, node.getYcoord()
                    + ( getTreeFontSet()._fm_large.getAscent() / down_shift_factor ), g );
        }
        if ( getControlPanel().isShowAnnotation() && node.getNodeData().isHasSequence()
                && ( node.getNodeData().getSequence().getAnnotations() != null )
                && ( !node.getNodeData().getSequence().getAnnotations().isEmpty() ) ) {
            if ( _sb.length() > 0 ) {
                x += getTreeFontSet()._fm_large.stringWidth( _sb.toString() ) + 5;
            }
            final Annotation ann = ( Annotation ) node.getNodeData().getSequence().getAnnotations().get( 0 );
            if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
                g.setColor( Color.BLACK );
            }
            else {
                g.setColor( calculateColorForAnnotation( ann ) );
            }
            final String ann_str = ann.asSimpleText().toString();
            TreePanel.drawString( ann_str, node.getXcoord() + x + 3 + TreePanel.HALF_BOX_SIZE, node.getYcoord()
                    + ( getTreeFontSet()._fm_large.getAscent() / down_shift_factor ), g );
            _sb.setLength( 0 );
            _sb.append( ann_str );
        }
        if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR )
                || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE )
                || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) ) {
            if ( ( getControlPanel().isShowBinaryCharacters() || getControlPanel().isShowBinaryCharacterCounts() )
                    && node.getNodeData().isHasBinaryCharacters() ) {
                if ( _sb.length() > 0 ) {
                    x += getTreeFontSet()._fm_large.stringWidth( _sb.toString() ) + 5;
                }
                if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
                    g.setColor( Color.BLACK );
                }
                else {
                    g.setColor( getTreeColorSet().getBinaryDomainCombinationsColor() );
                }
                if ( getControlPanel().isShowBinaryCharacters() ) {
                    TreePanel.drawString( node.getNodeData().getBinaryCharacters().getPresentCharactersAsStringBuffer()
                            .toString(), node.getXcoord() + x + 1 + TreePanel.HALF_BOX_SIZE, node.getYcoord()
                            + ( getTreeFontSet()._fm_large.getAscent() / down_shift_factor ), g );
                    paintGainedAndLostCharacters( g, node, node.getNodeData().getBinaryCharacters()
                            .getGainedCharactersAsStringBuffer().toString(), node.getNodeData().getBinaryCharacters()
                            .getLostCharactersAsStringBuffer().toString() );
                }
                else {
                    TreePanel.drawString( node.getNodeData().getBinaryCharacters().getPresentCount(), node.getXcoord()
                            + x + 2 + TreePanel.HALF_BOX_SIZE, node.getYcoord()
                            + ( getTreeFontSet()._fm_large.getAscent() / down_shift_factor ), g );
                    paintGainedAndLostCharacters( g, node, "+"
                            + node.getNodeData().getBinaryCharacters().getGainedCount(), "-"
                            + node.getNodeData().getBinaryCharacters().getLostCount() );
                }
            }
        }
    }

    final private void paintNodeDataUnrootedCirc( final Graphics2D g,
                                                  final PhylogenyNode node,
                                                  final boolean to_pdf,
                                                  final boolean to_graphics_file,
                                                  final boolean radial_labels,
                                                  final double ur_angle,
                                                  final boolean is_in_found_nodes ) {
        if ( isNodeDataInvisibleUnrootedCirc( node ) && !to_graphics_file && !to_pdf ) {
            return;
        }
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else if ( is_in_found_nodes ) {
            g.setColor( getTreeColorSet().getFoundColor() );
        }
        else if ( getControlPanel().isColorAccordingToTaxonomy() ) {
            g.setColor( getTaxonomyBasedColor( node ) );
        }
        else {
            g.setColor( getTreeColorSet().getSequenceColor() );
        }
        _sb.setLength( 0 );
        _sb.append( " " );
        if ( node.getNodeData().isHasTaxonomy()
                && ( getControlPanel().isShowTaxonomyCode() || getControlPanel().isShowTaxonomyNames() ) ) {
            final Taxonomy taxonomy = node.getNodeData().getTaxonomy();
            if ( _control_panel.isShowTaxonomyCode() && !ForesterUtil.isEmpty( taxonomy.getTaxonomyCode() ) ) {
                _sb.append( taxonomy.getTaxonomyCode() );
                _sb.append( " " );
            }
            if ( _control_panel.isShowTaxonomyNames() ) {
                if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() )
                        && !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                    _sb.append( taxonomy.getScientificName() );
                    _sb.append( " (" );
                    _sb.append( taxonomy.getCommonName() );
                    _sb.append( ") " );
                }
                else if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() ) ) {
                    _sb.append( taxonomy.getScientificName() );
                    _sb.append( " " );
                }
                else if ( !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                    _sb.append( taxonomy.getCommonName() );
                    _sb.append( " " );
                }
            }
        }
        if ( node.isCollapse() && ( ( !node.isRoot() && !node.getParent().isCollapse() ) || node.isRoot() ) ) {
            _sb.append( " [" );
            _sb.append( node.getAllExternalDescendants().size() );
            _sb.append( "]" );
        }
        if ( getControlPanel().isShowNodeNames() && ( node.getNodeName().length() > 0 ) ) {
            if ( _sb.length() > 0 ) {
                _sb.append( " " );
            }
            _sb.append( node.getNodeName() );
        }
        if ( node.getNodeData().isHasSequence() ) {
            if ( getControlPanel().isShowSequenceAcc() && ( node.getNodeData().getSequence().getAccession() != null ) ) {
                if ( _sb.length() > 0 ) {
                    _sb.append( " " );
                }
                if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() ) ) {
                    _sb.append( node.getNodeData().getSequence().getAccession().getSource() );
                    _sb.append( ":" );
                }
                _sb.append( node.getNodeData().getSequence().getAccession().getValue() );
            }
            if ( getControlPanel().isShowGeneNames() && ( node.getNodeData().getSequence().getName().length() > 0 ) ) {
                if ( _sb.length() > 0 ) {
                    _sb.append( " " );
                }
                _sb.append( node.getNodeData().getSequence().getName() );
            }
        }
        g.setFont( getTreeFontSet().getLargeFont() );
        if ( is_in_found_nodes ) {
            g.setFont( getTreeFontSet().getLargeFont().deriveFont( Font.BOLD ) );
        }
        if ( _sb.length() > 1 ) {
            final String sb_str = _sb.toString();
            double m = 0;
            if ( _graphics_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                m = _urt_nodeid_angle_map.get( node.getNodeId() ) % TWO_PI;
            }
            else {
                m = ( float ) ( ur_angle % TWO_PI );
            }
            _at = g.getTransform();
            boolean need_to_reset = false;
            final float x_coord = node.getXcoord();
            final float y_coord = node.getYcoord() + ( getTreeFontSet()._fm_large.getAscent() / 3.0f );
            if ( radial_labels ) {
                need_to_reset = true;
                boolean left = false;
                if ( ( m > HALF_PI ) && ( m < ONEHALF_PI ) ) {
                    m -= PI;
                    left = true;
                }
                g.rotate( m, x_coord, node.getYcoord() );
                if ( left ) {
                    g.translate( -( getTreeFontSet()._fm_large.getStringBounds( sb_str, g ).getWidth() ), 0 );
                }
            }
            else {
                if ( ( m > HALF_PI ) && ( m < ONEHALF_PI ) ) {
                    need_to_reset = true;
                    g.translate( -getTreeFontSet()._fm_large.getStringBounds( sb_str, g ).getWidth(), 0 );
                }
            }
            TreePanel.drawString( sb_str, x_coord, y_coord, g );
            if ( need_to_reset ) {
                g.setTransform( _at );
            }
        }
    }

    final private void paintNodeLite( final Graphics2D g, final PhylogenyNode node ) {
        if ( node.isCollapse() ) {
            if ( ( !node.isRoot() && !node.getParent().isCollapse() ) || node.isRoot() ) {
                paintCollapsedNode( g, node, false, false, false );
            }
            return;
        }
        if ( isInFoundNodes( node ) ) {
            g.setColor( getTreeColorSet().getFoundColor() );
            drawRectFilled( node.getXSecondary() - 1, node.getYSecondary() - 1, 3, 3, g );
        }
        float new_x = 0;
        if ( !node.isExternal() && !node.isCollapse() ) {
            boolean first_child = true;
            float y2 = 0.0f;
            final int parent_max_branch_to_leaf = getMaxBranchesToLeaf( node );
            for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
                final PhylogenyNode child_node = node.getChildNode( i );
                int factor_x;
                if ( !isUniformBranchLengthsForCladogram() ) {
                    factor_x = node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes();
                }
                else {
                    factor_x = parent_max_branch_to_leaf - getMaxBranchesToLeaf( child_node );
                }
                if ( first_child ) {
                    first_child = false;
                    y2 = node.getYSecondary()
                            - ( getOvYDistance() * ( node.getNumberOfExternalNodes() - child_node
                                    .getNumberOfExternalNodes() ) );
                }
                else {
                    y2 += getOvYDistance() * child_node.getNumberOfExternalNodes();
                }
                final float x2 = calculateOvBranchLengthToParent( child_node, factor_x );
                new_x = x2 + node.getXSecondary();
                final float diff_y = node.getYSecondary() - y2;
                final float diff_x = node.getXSecondary() - new_x;
                if ( ( diff_y > 2 ) || ( diff_y < -2 ) || ( diff_x > 2 ) || ( diff_x < -2 ) ) {
                    paintBranchLite( g, node.getXSecondary(), new_x, node.getYSecondary(), y2, child_node );
                }
                child_node.setXSecondary( new_x );
                child_node.setYSecondary( y2 );
                y2 += getOvYDistance() * child_node.getNumberOfExternalNodes();
            }
        }
    }

    final private void paintNodeRectangular( final Graphics2D g,
                                             final PhylogenyNode node,
                                             final boolean to_pdf,
                                             final boolean dynamically_hide,
                                             final int dynamic_hiding_factor,
                                             final boolean to_graphics_file ) {
        final boolean is_in_found_nodes = isInFoundNodes( node );
        if ( node.isCollapse() ) {
            if ( ( !node.isRoot() && !node.getParent().isCollapse() ) || node.isRoot() ) {
                paintCollapsedNode( g, node, to_graphics_file, to_pdf, is_in_found_nodes );
            }
            return;
        }
        if ( node.isExternal() ) {
            ++_external_node_index;
        }
        // Confidence values
        if ( getControlPanel().isShowBootstrapValues()
                && !node.isExternal()
                && !node.isRoot()
                && ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED )
                        || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) )
                && node.getBranchData().isHasConfidences() ) {
            paintConfidenceValues( g, node, to_pdf, to_graphics_file );
        }
        // Draw a line to root:
        if ( node.isRoot() && _phylogeny.isRooted() ) {
            paintRootBranch( g, node.getXcoord(), node.getYcoord(), node, to_pdf, to_graphics_file );
        }
        float new_x = 0;
        float new_x_min = Float.MAX_VALUE;
        final boolean disallow_shortcutting = dynamic_hiding_factor < 40;
        float min_dist = 1.5f;
        if ( !disallow_shortcutting ) {
            //   System.out.println( dynamic_hiding_factor );
            if ( dynamic_hiding_factor > 4000 ) {
                min_dist = 4;
            }
            else if ( dynamic_hiding_factor > 1000 ) {
                min_dist = 3;
            }
            else if ( dynamic_hiding_factor > 100 ) {
                min_dist = 2;
            }
        }
        if ( !node.isExternal() && !node.isCollapse() ) {
            boolean first_child = true;
            float y2 = 0.0f;
            final int parent_max_branch_to_leaf = getMaxBranchesToLeaf( node );
            for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
                final PhylogenyNode child_node = node.getChildNode( i );
                int factor_x;
                if ( !isUniformBranchLengthsForCladogram() ) {
                    factor_x = node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes();
                }
                else {
                    factor_x = parent_max_branch_to_leaf - getMaxBranchesToLeaf( child_node );
                }
                if ( first_child ) {
                    first_child = false;
                    y2 = node.getYcoord()
                            - ( _y_distance * ( node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes() ) );
                }
                else {
                    y2 += _y_distance * child_node.getNumberOfExternalNodes();
                }
                final float x2 = calculateBranchLengthToParent( child_node, factor_x );
                new_x = x2 + node.getXcoord();
                if ( dynamically_hide && ( x2 < new_x_min ) ) {
                    new_x_min = x2;
                }
                final float diff_y = node.getYcoord() - y2;
                final float diff_x = node.getXcoord() - new_x;
                if ( disallow_shortcutting || ( diff_y > min_dist ) || ( diff_y < -min_dist ) || ( diff_x > min_dist )
                        || ( diff_x < -min_dist ) || to_graphics_file || to_pdf ) {
                    paintBranchRectangular( g,
                                            node.getXcoord(),
                                            new_x,
                                            node.getYcoord(),
                                            y2,
                                            child_node,
                                            to_pdf,
                                            to_graphics_file );
                }
                child_node.setXcoord( new_x );
                child_node.setYcoord( y2 );
                y2 += _y_distance * child_node.getNumberOfExternalNodes();
            }
        }
        if ( dynamically_hide
                && !is_in_found_nodes
                && ( ( node.isExternal() && ( _external_node_index % dynamic_hiding_factor != 1 ) ) || ( !node
                        .isExternal() && ( ( new_x_min < 20 ) || ( _y_distance * node.getNumberOfExternalNodes() < getTreeFontSet()._fm_large
                        .getHeight() ) ) ) ) ) {
            return;
        }
        paintNodeData( g, node, to_graphics_file, to_pdf, is_in_found_nodes );
        paintNodeWithRenderableData( g, node, to_graphics_file, to_pdf );
    }

    final private void paintNodeWithRenderableData( final Graphics2D g,
                                                    final PhylogenyNode node,
                                                    final boolean to_graphics_file,
                                                    final boolean to_pdf ) {
        if ( isNodeDataInvisible( node ) && !to_graphics_file ) {
            return;
        }
        if ( ( !getControlPanel().isShowInternalData() && !node.isExternal() ) ) {
            return;
        }
        if ( getControlPanel().isShowDomainArchitectures() && node.getNodeData().isHasSequence()
                && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
            RenderableDomainArchitecture rds = null;
            try {
                rds = ( RenderableDomainArchitecture ) node.getNodeData().getSequence().getDomainArchitecture();
            }
            catch ( final ClassCastException cce ) {
                return;
            }
            rds.setRenderingHeight( 6 );
            int x = 0;
            if ( getControlPanel().isShowTaxonomyCode() && ( PhylogenyMethods.getSpecies( node ).length() > 0 ) ) {
                x += getTreeFontSet()._fm_large_italic.stringWidth( PhylogenyMethods.getSpecies( node ) + " " );
            }
            if ( getControlPanel().isShowNodeNames() && ( node.getNodeName().length() > 0 ) ) {
                x += getTreeFontSet()._fm_large.stringWidth( node.getNodeName() + " " );
            }
            rds.render( node.getXcoord() + x, node.getYcoord() - 3, g, this, to_pdf );
        }
    }

    final private void paintOvRectangle( final Graphics2D g ) {
        final float w_ratio = ( float ) getWidth() / getVisibleRect().width;
        final float h_ratio = ( float ) getHeight() / getVisibleRect().height;
        final float x_ratio = ( float ) getWidth() / getVisibleRect().x;
        final float y_ratio = ( float ) getHeight() / getVisibleRect().y;
        final float width = getOvMaxWidth() / w_ratio;
        final float height = getOvMaxHeight() / h_ratio;
        final float x = getVisibleRect().x + getOvXPosition() + getOvMaxWidth() / x_ratio;
        final float y = getVisibleRect().y + getOvYPosition() + getOvMaxHeight() / y_ratio;
        g.setColor( getTreeColorSet().getFoundColor() );
        getOvRectangle().setRect( x, y, width, height );
        if ( ( width < 6 ) && ( height < 6 ) ) {
            drawRectFilled( x, y, 6, 6, g );
            getOvVirtualRectangle().setRect( x, y, 6, 6 );
        }
        else if ( width < 6 ) {
            drawRectFilled( x, y, 6, height, g );
            getOvVirtualRectangle().setRect( x, y, 6, height );
        }
        else if ( height < 6 ) {
            drawRectFilled( x, y, width, 6, g );
            getOvVirtualRectangle().setRect( x, y, width, 6 );
        }
        else {
            drawRect( x, y, width, height, g );
            if ( isInOvRect() ) {
                drawRect( x + 1, y + 1, width - 2, height - 2, g );
            }
            getOvVirtualRectangle().setRect( x, y, width, height );
        }
    }

    final void paintPhylogeny( final Graphics2D g,
                               final boolean to_pdf,
                               final boolean to_graphics_file,
                               final int graphics_file_width,
                               final int graphics_file_height,
                               final int graphics_file_x,
                               final int graphics_file_y ) {
        // Color the background
        if ( !to_pdf ) {
            final Rectangle r = getVisibleRect();
            if ( !getOptions().isBackgroundColorGradient() || getOptions().isPrintBlackAndWhite() ) {
                g.setColor( getTreeColorSet().getBackgroundColor() );
                if ( !to_graphics_file ) {
                    g.fill( r );
                }
                else {
                    if ( getOptions().isPrintBlackAndWhite() ) {
                        g.setColor( Color.WHITE );
                    }
                    g.fillRect( graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height );
                }
            }
            else {
                if ( !to_graphics_file ) {
                    g.setPaint( new GradientPaint( r.x, r.y, getTreeColorSet().getBackgroundColor(), r.x, r.y
                            + r.height, getTreeColorSet().getBackgroundColorGradientBottom() ) );
                    g.fill( r );
                }
                else {
                    g.setPaint( new GradientPaint( graphics_file_x,
                                                   graphics_file_y,
                                                   getTreeColorSet().getBackgroundColor(),
                                                   graphics_file_x,
                                                   graphics_file_y + graphics_file_height,
                                                   getTreeColorSet().getBackgroundColorGradientBottom() ) );
                    g.fillRect( graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height );
                }
            }
            g.setStroke( new BasicStroke( 1 ) );
        }
        else {
            g.setStroke( new BasicStroke( getOptions().getPrintLineWidth() ) );
        }
        if ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                && ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
            _external_node_index = 0;
            // Position starting X of tree
            if ( !_phylogeny.isRooted() ) {
                _phylogeny.getRoot().setXcoord( TreePanel.MOVE );
            }
            else if ( ( _phylogeny.getRoot().getDistanceToParent() > 0.0 ) && getControlPanel().isDrawPhylogram() ) {
                _phylogeny.getRoot().setXcoord( ( float ) ( TreePanel.MOVE + ( _phylogeny.getRoot()
                        .getDistanceToParent() * getXcorrectionFactor() ) ) );
            }
            else {
                _phylogeny.getRoot().setXcoord( TreePanel.MOVE + getXdistance() );
            }
            // Position starting Y of tree
            _phylogeny.getRoot().setYcoord( ( getYdistance() * _phylogeny.getRoot().getNumberOfExternalNodes() )
                    + ( TreePanel.MOVE / 2.0f ) );
            final int dynamic_hiding_factor = ( int ) ( getTreeFontSet()._fm_large.getHeight() / ( 1.5 * getYdistance() ) );
            if ( getControlPanel().isDynamicallyHideData() ) {
                if ( dynamic_hiding_factor > 1 ) {
                    getControlPanel().setDynamicHidingIsOn( true );
                }
                else {
                    getControlPanel().setDynamicHidingIsOn( false );
                }
            }
            final PhylogenyNodeIterator it;
            for( it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
                paintNodeRectangular( g, it.next(), to_pdf, getControlPanel().isDynamicallyHideData()
                        && ( dynamic_hiding_factor > 1 ), dynamic_hiding_factor, to_graphics_file );
            }
            if ( getOptions().isShowScale() ) {
                if ( !( to_graphics_file || to_pdf ) ) {
                    paintScale( g,
                                getVisibleRect().x,
                                getVisibleRect().y + getVisibleRect().height,
                                to_pdf,
                                to_graphics_file );
                }
                else {
                    paintScale( g, graphics_file_x, graphics_file_y + graphics_file_height, to_pdf, to_graphics_file );
                }
            }
            if ( getOptions().isShowOverview() && isOvOn() && !to_graphics_file && !to_pdf ) {
                paintPhylogenyLite( g );
            }
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            if ( getControlPanel().getDynamicallyHideData() != null ) {
                getControlPanel().setDynamicHidingIsOn( false );
            }
            final double angle = getStartingAngle();
            final boolean radial_labels = getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL;
            _dynamic_hiding_factor = 0;
            if ( getControlPanel().isDynamicallyHideData() ) {
                _dynamic_hiding_factor = ( int ) ( ( getTreeFontSet()._fm_large.getHeight() * 1.5 * getPhylogeny()
                        .getNumberOfExternalNodes() ) / ( TWO_PI * 10 ) );
            }
            if ( getControlPanel().getDynamicallyHideData() != null ) {
                if ( _dynamic_hiding_factor > 1 ) {
                    getControlPanel().setDynamicHidingIsOn( true );
                }
                else {
                    getControlPanel().setDynamicHidingIsOn( false );
                }
            }
            paintUnrooted( _phylogeny.getRoot(),
                           angle,
                           ( float ) ( angle + 2 * Math.PI ),
                           radial_labels,
                           g,
                           to_pdf,
                           to_graphics_file );
            if ( getOptions().isShowScale() ) {
                if ( !( to_graphics_file || to_pdf ) ) {
                    paintScale( g,
                                getVisibleRect().x,
                                getVisibleRect().y + getVisibleRect().height,
                                to_pdf,
                                to_graphics_file );
                }
                else {
                    paintScale( g, graphics_file_x, graphics_file_y + graphics_file_height, to_pdf, to_graphics_file );
                }
            }
            if ( getOptions().isShowOverview() && isOvOn() && !to_graphics_file && !to_pdf ) {
                g.setColor( getTreeColorSet().getOvColor() );
                paintUnrootedLite( _phylogeny.getRoot(),
                                   angle,
                                   angle + 2 * Math.PI,
                                   g,
                                   ( getUrtFactorOv() / ( getVisibleRect().width / getOvMaxWidth() ) ) );
                paintOvRectangle( g );
            }
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
            final int radius = ( int ) ( ( Math.min( getPreferredSize().getWidth(), getPreferredSize().getHeight() ) / 2 ) - ( MOVE + getLongestExtNodeInfo() ) );
            final int d = radius + MOVE + getLongestExtNodeInfo();
            _dynamic_hiding_factor = 0;
            if ( getControlPanel().isDynamicallyHideData() && ( radius > 0 ) ) {
                _dynamic_hiding_factor = ( int ) ( ( getTreeFontSet()._fm_large.getHeight() * 1.5 * getPhylogeny()
                        .getNumberOfExternalNodes() ) / ( TWO_PI * radius ) );
            }
            if ( getControlPanel().getDynamicallyHideData() != null ) {
                if ( _dynamic_hiding_factor > 1 ) {
                    getControlPanel().setDynamicHidingIsOn( true );
                }
                else {
                    getControlPanel().setDynamicHidingIsOn( false );
                }
            }
            paintCircular( _phylogeny, getStartingAngle(), d, d, radius > 0 ? radius : 0, g, to_pdf, to_graphics_file );
            if ( getOptions().isShowOverview() && isOvOn() && !to_graphics_file && !to_pdf ) {
                final int radius_ov = ( int ) ( getOvMaxHeight() < getOvMaxWidth() ? getOvMaxHeight() / 2
                        : getOvMaxWidth() / 2 );
                double x_scale = 1.0;
                double y_scale = 1.0;
                int x_pos = getVisibleRect().x + getOvXPosition();
                int y_pos = getVisibleRect().y + getOvYPosition();
                if ( getWidth() > getHeight() ) {
                    x_scale = ( double ) getHeight() / getWidth();
                    x_pos = ForesterUtil.roundToInt( x_pos / x_scale );
                }
                else {
                    y_scale = ( double ) getWidth() / getHeight();
                    y_pos = ForesterUtil.roundToInt( y_pos / y_scale );
                }
                _at = g.getTransform();
                g.scale( x_scale, y_scale );
                paintCircularLite( _phylogeny,
                                   getStartingAngle(),
                                   x_pos + radius_ov,
                                   y_pos + radius_ov,
                                   ( int ) ( radius_ov - ( getLongestExtNodeInfo() / ( getVisibleRect().width / getOvRectangle()
                                           .getWidth() ) ) ),
                                   g );
                g.setTransform( _at );
                paintOvRectangle( g );
            }
        }
    }

    final private void paintPhylogenyLite( final Graphics2D g ) {
        _phylogeny
                .getRoot()
                .setXSecondary( ( float ) ( getVisibleRect().x + getOvXPosition() + ( MOVE / ( getVisibleRect().width / getOvRectangle()
                        .getWidth() ) ) ) );
        _phylogeny.getRoot().setYSecondary( ( getVisibleRect().y + getOvYStart() ) );
        final PhylogenyNodeIterator it;
        for( it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
            paintNodeLite( g, it.next() );
        }
        paintOvRectangle( g );
    }

    /**
     * Paint the root branch. (Differs from others because it will always be a
     * single horizontal line).
     * @param to_graphics_file 
     * 
     * @return new x1 value
     */
    final private void paintRootBranch( final Graphics2D g,
                                        final float x1,
                                        final float y1,
                                        final PhylogenyNode root,
                                        final boolean to_pdf,
                                        final boolean to_graphics_file ) {
        assignGraphicsForBranchWithColorForParentBranch( root, false, g, to_pdf, to_graphics_file );
        float d = getXdistance();
        if ( getControlPanel().isDrawPhylogram() && ( root.getDistanceToParent() > 0.0 ) ) {
            d = ( float ) ( getXcorrectionFactor() * root.getDistanceToParent() );
        }
        if ( d < MIN_ROOT_LENGTH ) {
            d = MIN_ROOT_LENGTH;
        }
        if ( !getControlPanel().isWidthBranches() || ( PhylogenyMethods.getBranchWidthValue( root ) == 1 ) ) {
            drawLine( x1 - d, root.getYcoord(), x1, root.getYcoord(), g );
        }
        else {
            final double w = PhylogenyMethods.getBranchWidthValue( root );
            drawRectFilled( x1 - d, root.getYcoord() - ( w / 2 ), d, w, g );
        }
        paintNodeBox( x1, root.getYcoord(), root, g, to_pdf, to_graphics_file, isInFoundNodes( root ) );
    }

    final private void paintScale( final Graphics2D g,
                                   int x1,
                                   int y1,
                                   final boolean to_pdf,
                                   final boolean to_graphics_file ) {
        if ( !getControlPanel().isDrawPhylogram() || ( getScaleDistance() <= 0.0 ) ) {
            return;
        }
        x1 += MOVE;
        final double x2 = x1 + ( getScaleDistance() * getXcorrectionFactor() );
        y1 -= 12;
        final int y2 = y1 - 8;
        final int y3 = y1 - 4;
        g.setFont( getTreeFontSet().getSmallFont() );
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else {
            g.setColor( getTreeColorSet().getBranchLengthColor() );
        }
        drawLine( x1, y1, x1, y2, g );
        drawLine( x2, y1, x2, y2, g );
        drawLine( x1, y3, x2, y3, g );
        if ( getScaleLabel() != null ) {
            g.drawString( getScaleLabel(), ( x1 + 2 ), y3 - 2 );
        }
    }

    final private int paintTaxonomy( final Graphics2D g,
                                     final PhylogenyNode node,
                                     final boolean is_in_found_nodes,
                                     final boolean to_pdf,
                                     final boolean to_graphics_file ) {
        final Taxonomy taxonomy = node.getNodeData().getTaxonomy();
        g.setFont( getTreeFontSet().getLargeItalicFont() );
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else if ( is_in_found_nodes ) {
            g.setFont( getTreeFontSet().getLargeItalicFont().deriveFont( TreeFontSet.BOLD_AND_ITALIC ) );
            g.setColor( getTreeColorSet().getFoundColor() );
        }
        else if ( getControlPanel().isColorAccordingToTaxonomy() ) {
            g.setColor( getTaxonomyBasedColor( node ) );
        }
        else {
            g.setColor( getTreeColorSet().getTaxonomyColor() );
        }
        final double start_x = node.getXcoord() + 3 + TreePanel.HALF_BOX_SIZE;
        final double start_y = node.getYcoord()
                + ( getTreeFontSet()._fm_large.getAscent() / ( node.getNumberOfDescendants() == 1 ? 1 : 3.0 ) );
        _sb.setLength( 0 );
        if ( _control_panel.isShowTaxonomyCode() && !ForesterUtil.isEmpty( taxonomy.getTaxonomyCode() ) ) {
            _sb.append( taxonomy.getTaxonomyCode() );
            _sb.append( " " );
        }
        if ( _control_panel.isShowTaxonomyNames() ) {
            if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() )
                    && !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                _sb.append( taxonomy.getScientificName() );
                _sb.append( " (" );
                _sb.append( taxonomy.getCommonName() );
                _sb.append( ") " );
            }
            else if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() ) ) {
                _sb.append( taxonomy.getScientificName() );
                _sb.append( " " );
            }
            else if ( !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                _sb.append( taxonomy.getCommonName() );
                _sb.append( " " );
            }
        }
        final String label = _sb.toString();
        TreePanel.drawString( label, start_x, start_y, g );
        if ( is_in_found_nodes ) {
            return getTreeFontSet()._fm_large_italic_bold.stringWidth( label );
        }
        else {
            return getTreeFontSet()._fm_large_italic.stringWidth( label );
        }
    }

    final private void paintUnrooted( final PhylogenyNode n,
                                      final double low_angle,
                                      final double high_angle,
                                      final boolean radial_labels,
                                      final Graphics2D g,
                                      final boolean to_pdf,
                                      final boolean to_graphics_file ) {
        if ( n.isRoot() ) {
            n.setXcoord( getWidth() / 2 );
            n.setYcoord( getHeight() / 2 );
            paintNodeBox( n.getXcoord(), n.getYcoord(), n, g, to_pdf, to_graphics_file, isInFoundNodes( n ) );
        }
        if ( n.isExternal() ) {
            paintNodeDataUnrootedCirc( g,
                                       n,
                                       to_pdf,
                                       to_graphics_file,
                                       radial_labels,
                                       ( high_angle + low_angle ) / 2,
                                       isInFoundNodes( n ) );
            return;
        }
        final float num_enclosed = n.getNumberOfExternalNodes();
        final float x = n.getXcoord();
        final float y = n.getYcoord();
        double current_angle = low_angle;
        // final boolean n_below = n.getYcoord() < getVisibleRect().getMinY() - 20;
        // final boolean n_above = n.getYcoord() > getVisibleRect().getMaxY() + 20;
        // final boolean n_left = n.getXcoord() < getVisibleRect().getMinX() - 20;
        // final boolean n_right = n.getXcoord() > getVisibleRect().getMaxX() + 20;
        for( int i = 0; i < n.getNumberOfDescendants(); ++i ) {
            final PhylogenyNode desc = n.getChildNode( i );
            ///  if ( ( ( n_below ) & ( desc.getYcoord() < getVisibleRect().getMinY() - 20 ) )
            //          || ( ( n_above ) & ( desc.getYcoord() > getVisibleRect().getMaxY() + 20 ) )
            //         || ( ( n_left ) & ( desc.getXcoord() < getVisibleRect().getMinX() - 20 ) )
            //          || ( ( n_right ) & ( desc.getXcoord() > getVisibleRect().getMaxX() + 20 ) ) ) {
            //     continue;
            // }
            //if ( ( desc.getYcoord() > n.getYcoord() ) && ( n.getYcoord() > getVisibleRect().getMaxY() - 20 ) ) {
            //    continue;
            //}
            //if ( ( desc.getYcoord() < n.getYcoord() ) && ( n.getYcoord() < getVisibleRect().getMinY() + 20 ) ) {
            //    continue;
            // }
            final int desc_num_enclosed = desc.getNumberOfExternalNodes();
            final double arc_size = ( desc_num_enclosed / num_enclosed ) * ( high_angle - low_angle );
            float length;
            if ( isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram() ) {
                if ( desc.getDistanceToParent() < 0 ) {
                    length = 0;
                }
                else {
                    length = ( float ) ( desc.getDistanceToParent() * getUrtFactor() );
                }
            }
            else {
                length = getUrtFactor();
            }
            final double mid_angle = current_angle + arc_size / 2;
            final float new_x = ( float ) ( x + Math.cos( mid_angle ) * length );
            final float new_y = ( float ) ( y + Math.sin( mid_angle ) * length );
            desc.setXcoord( new_x );
            desc.setYcoord( new_y );
            paintNodeBox( new_x, new_y, desc, g, to_pdf, to_graphics_file, isInFoundNodes( desc ) );
            paintUnrooted( desc, current_angle, current_angle + arc_size, radial_labels, g, to_pdf, to_graphics_file );
            current_angle += arc_size;
            assignGraphicsForBranchWithColorForParentBranch( desc, false, g, to_pdf, to_graphics_file );
            drawLine( x, y, new_x, new_y, g );
        }
    }

    final private void paintUnrootedLite( final PhylogenyNode n,
                                          final double low_angle,
                                          final double high_angle,
                                          final Graphics2D g,
                                          final float urt_ov_factor ) {
        if ( n.isRoot() ) {
            final int x_pos = ( int ) ( getVisibleRect().x + getOvXPosition() + getOvMaxWidth() / 2 );
            final int y_pos = ( int ) ( getVisibleRect().y + getOvYPosition() + getOvMaxHeight() / 2 );
            n.setXSecondary( x_pos );
            n.setYSecondary( y_pos );
        }
        if ( n.isExternal() ) {
            return;
        }
        final float num_enclosed = n.getNumberOfExternalNodes();
        final float x = n.getXSecondary();
        final float y = n.getYSecondary();
        double current_angle = low_angle;
        for( int i = 0; i < n.getNumberOfDescendants(); ++i ) {
            final PhylogenyNode desc = n.getChildNode( i );
            final int desc_num_enclosed = desc.getNumberOfExternalNodes();
            final double arc_size = ( desc_num_enclosed / num_enclosed ) * ( high_angle - low_angle );
            float length;
            if ( isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram() ) {
                if ( desc.getDistanceToParent() < 0 ) {
                    length = 0;
                }
                else {
                    length = ( float ) ( desc.getDistanceToParent() * urt_ov_factor );
                }
            }
            else {
                length = urt_ov_factor;
            }
            final double mid_angle = current_angle + arc_size / 2;
            final float new_x = ( float ) ( x + Math.cos( mid_angle ) * length );
            final float new_y = ( float ) ( y + Math.sin( mid_angle ) * length );
            desc.setXSecondary( new_x );
            desc.setYSecondary( new_y );
            if ( isInFoundNodes( desc ) ) {
                g.setColor( getTreeColorSet().getFoundColor() );
                drawRectFilled( desc.getXSecondary() - 1, desc.getYSecondary() - 1, 3, 3, g );
                g.setColor( getTreeColorSet().getOvColor() );
            }
            paintUnrootedLite( desc, current_angle, current_angle + arc_size, g, urt_ov_factor );
            current_angle += arc_size;
            drawLine( x, y, new_x, new_y, g );
        }
    }

    final private void pasteSubtree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        if ( ( getCutOrCopiedTree() == null ) || getCutOrCopiedTree().isEmpty() ) {
            JOptionPane.showMessageDialog( this,
                                           "No tree in buffer (need to copy or cut a subtree first)",
                                           "Attempt to paste with empty buffer",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        final String label = getASimpleTextRepresentationOfANode( getCutOrCopiedTree().getRoot() );
        final Object[] options = { "As sibling", "As descendant", "Cancel" };
        final int r = JOptionPane.showOptionDialog( this,
                                                    "How to paste subtree" + label + "?",
                                                    "Paste Subtree",
                                                    JOptionPane.CLOSED_OPTION,
                                                    JOptionPane.QUESTION_MESSAGE,
                                                    null,
                                                    options,
                                                    options[ 2 ] );
        boolean paste_as_sibling = true;
        if ( r == 1 ) {
            paste_as_sibling = false;
        }
        else if ( r != 0 ) {
            return;
        }
        final Phylogeny buffer_phy = getCutOrCopiedTree().copy();
        buffer_phy.setAllNodesToNotCollapse();
        buffer_phy.preOrderReId();
        if ( paste_as_sibling ) {
            if ( node.isRoot() ) {
                JOptionPane.showMessageDialog( this,
                                               "Cannot paste sibling to root",
                                               "Attempt to paste sibling to root",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            buffer_phy.addAsSibling( node );
        }
        else {
            buffer_phy.addAsChild( node );
        }
        if ( getCopiedAndPastedNodes() == null ) {
            setCopiedAndPastedNodes( new HashSet<PhylogenyNode>() );
        }
        getCopiedAndPastedNodes().addAll( PhylogenyMethods.obtainAllNodesAsSet( buffer_phy ) );
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.hashIDs();
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        setEdited( true );
        repaint();
    }

    final public int print( final Graphics g, final PageFormat page_format, final int page_index )
            throws PrinterException {
        if ( page_index > 0 ) {
            return ( NO_SUCH_PAGE );
        }
        else {
            final Graphics2D g2d = ( Graphics2D ) g;
            g2d.translate( page_format.getImageableX(), page_format.getImageableY() );
            // Turn off double buffering !?
            paintPhylogeny( g2d, true, false, 0, 0, 0, 0 );
            // Turn double buffering back on !?
            return ( PAGE_EXISTS );
        }
    }

    final void recalculateMaxDistanceToRoot() {
        _max_distance_to_root = PhylogenyMethods.calculateMaxDistanceToRoot( getPhylogeny() );
    }

    /**
     * Remove all edit-node frames
     */
    final void removeAllEditNodeJFrames() {
        for( int i = 0; i <= ( TreePanel.MAX_NODE_FRAMES - 1 ); i++ ) {
            if ( _node_frames[ i ] != null ) {
                _node_frames[ i ].dispose();
                _node_frames[ i ] = null;
            }
        }
        _node_frame_index = 0;
    }

    /**
     * Remove a node-edit frame.
     */
    final void removeEditNodeFrame( final int i ) {
        _node_frame_index--;
        _node_frames[ i ] = null;
        if ( i < _node_frame_index ) {
            for( int j = 0; j < _node_frame_index - 1; j++ ) {
                _node_frames[ j ] = _node_frames[ j + 1 ];
            }
            _node_frames[ _node_frame_index ] = null;
        }
    }

    final void reRoot( final PhylogenyNode node ) {
        if ( !getPhylogeny().isRerootable() ) {
            JOptionPane.showMessageDialog( this,
                                           "This is not rerootable",
                                           "Not rerootable",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot reroot in unrooted display type",
                                           "Attempt to reroot tree in unrooted display",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        getPhylogeny().reRoot( node );
        getPhylogeny().recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        resetPreferredSize();
        getMainPanel().adjustJScrollPane();
        repaint();
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
            getControlPanel().showWhole();
        }
    }

    final void resetNodeIdToDistToLeafMap() {
        _nodeid_dist_to_leaf = new HashMap<Integer, Short>();
    }

    final void resetPreferredSize() {
        if ( ( getPhylogeny() == null ) || getPhylogeny().isEmpty() ) {
            return;
        }
        int x = 0;
        int y = 0;
        y = TreePanel.MOVE
                + ForesterUtil.roundToInt( getYdistance() * getPhylogeny().getRoot().getNumberOfExternalNodes() * 2 );
        if ( getControlPanel().isDrawPhylogram() ) {
            x = TreePanel.MOVE
                    + getLongestExtNodeInfo()
                    + ForesterUtil
                            .roundToInt( ( getXcorrectionFactor() * getPhylogeny().getHeight() ) + getXdistance() );
        }
        else {
            if ( !isNonLinedUpCladogram() && !isUniformBranchLengthsForCladogram() ) {
                x = TreePanel.MOVE
                        + getLongestExtNodeInfo()
                        + ForesterUtil.roundToInt( getXdistance()
                                * ( getPhylogeny().getRoot().getNumberOfExternalNodes() + 2 ) );
            }
            else {
                x = TreePanel.MOVE
                        + getLongestExtNodeInfo()
                        + ForesterUtil.roundToInt( getXdistance()
                                * ( PhylogenyMethods.calculateMaxDepth( getPhylogeny() ) + 1 ) );
            }
        }
        setPreferredSize( new Dimension( x, y ) );
    }

    final void setArrowCursor() {
        setCursor( ARROW_CURSOR );
        repaint();
    }

    final void setControlPanel( final ControlPanel atv_control ) {
        _control_panel = atv_control;
    }

    final private void setCopiedAndPastedNodes( final Set<PhylogenyNode> copied_and_pasted_nodes ) {
        getMainPanel().setCopiedAndPastedNodes( copied_and_pasted_nodes );
    }

    final private void setCutOrCopiedTree( final Phylogeny cut_or_copied_tree ) {
        getMainPanel().setCutOrCopiedTree( cut_or_copied_tree );
    }

    final void setEdited( final boolean edited ) {
        _edited = edited;
    }

    final void setFoundNodes( final Set<PhylogenyNode> found_nodes ) {
        _found_nodes = found_nodes;
    }

    final private void setInOv( final boolean in_ov ) {
        _in_ov = in_ov;
    }

    final void setInOvRect( final boolean in_ov_rect ) {
        _in_ov_rect = in_ov_rect;
    }

    final void setLargeFonts() {
        getTreeFontSet().largeFonts();
    }

    final void setLastMouseDragPointX( final float x ) {
        _last_drag_point_x = x;
    }

    final void setLastMouseDragPointY( final float y ) {
        _last_drag_point_y = y;
    }

    final void setLongestExtNodeInfo( final int i ) {
        _longest_ext_node_info = i;
    }

    final void setMediumFonts() {
        getTreeFontSet().mediumFonts();
    }

    final private void setOvMaxHeight( final float ov_max_height ) {
        _ov_max_height = ov_max_height;
    }

    final private void setOvMaxWidth( final float ov_max_width ) {
        _ov_max_width = ov_max_width;
    }

    final void setOvOn( final boolean ov_on ) {
        _ov_on = ov_on;
    }

    final private void setOvXcorrectionFactor( final float f ) {
        _ov_x_correction_factor = f;
    }

    final private void setOvXDistance( final float ov_x_distance ) {
        _ov_x_distance = ov_x_distance;
    }

    final private void setOvXPosition( final int ov_x_position ) {
        _ov_x_position = ov_x_position;
    }

    final private void setOvYDistance( final float ov_y_distance ) {
        _ov_y_distance = ov_y_distance;
    }

    final private void setOvYPosition( final int ov_y_position ) {
        _ov_y_position = ov_y_position;
    }

    final private void setOvYStart( final int ov_y_start ) {
        _ov_y_start = ov_y_start;
    }

    /**
     * Set parameters for printing the displayed tree
     * 
     * @param x
     * @param y
     */
    final void setParametersForPainting( final int x, final int y, final boolean recalc_longest_ext_node_info ) {
        // updateStyle(); not needed?
        if ( ( _phylogeny != null ) && !_phylogeny.isEmpty() ) {
            initNodeData();
            if ( recalc_longest_ext_node_info ) {
                calculateLongestExtNodeInfo();
            }
            int ext_nodes = _phylogeny.getRoot().getNumberOfExternalNodes();
            final int max_depth = PhylogenyMethods.calculateMaxDepth( _phylogeny );
            if ( ext_nodes == 1 ) {
                ext_nodes = max_depth;
                if ( ext_nodes < 1 ) {
                    ext_nodes = 1;
                }
            }
            updateOvSizes();
            float xdist = 0;
            float ov_xdist = 0;
            if ( !isNonLinedUpCladogram() && !isUniformBranchLengthsForCladogram() ) {
                xdist = ( float ) ( ( x - getLongestExtNodeInfo() - TreePanel.MOVE ) / ( ext_nodes + 3.0 ) );
                ov_xdist = ( float ) ( getOvMaxWidth() / ( ext_nodes + 3.0 ) );
            }
            else {
                xdist = ( ( x - getLongestExtNodeInfo() - TreePanel.MOVE ) / ( max_depth + 1 ) );
                ov_xdist = ( getOvMaxWidth() / ( max_depth + 1 ) );
            }
            float ydist = ( float ) ( ( y - TreePanel.MOVE ) / ( ext_nodes * 2.0 ) );
            if ( xdist < 0.0 ) {
                xdist = 0.0f;
            }
            if ( ov_xdist < 0.0 ) {
                ov_xdist = 0.0f;
            }
            if ( ydist < 0.0 ) {
                ydist = 0.0f;
            }
            setXdistance( xdist );
            setYdistance( ydist );
            setOvXDistance( ov_xdist );
            final double height = _phylogeny.getHeight();
            if ( height > 0 ) {
                final float corr = ( float ) ( ( x - TreePanel.MOVE - getLongestExtNodeInfo() - getXdistance() ) / height );
                setXcorrectionFactor( corr > 0 ? corr : 0 );
                final float ov_corr = ( float ) ( ( getOvMaxWidth() - getOvXDistance() ) / height );
                setOvXcorrectionFactor( ov_corr > 0 ? ov_corr : 0 );
            }
            else {
                setXcorrectionFactor( 0 );
                setOvXcorrectionFactor( 0 );
            }
            _circ_max_depth = max_depth;
            setUpUrtFactor();
        }
    }

    final void setPhylogenyGraphicsType( final PHYLOGENY_GRAPHICS_TYPE graphics_type ) {
        _graphics_type = graphics_type;
        setTextAntialias();
    }

    final private void setScaleDistance( final double scale_distance ) {
        _scale_distance = scale_distance;
    }

    final private void setScaleLabel( final String scale_label ) {
        _scale_label = scale_label;
    }

    final void setSmallFonts() {
        getTreeFontSet().smallFonts();
    }

    final void setStartingAngle( final double starting_angle ) {
        _urt_starting_angle = starting_angle;
    }

    final void setSuperTinyFonts() {
        getTreeFontSet().superTinyFonts();
    }

    final void setTextAntialias() {
        if ( ( _phylogeny != null ) && !_phylogeny.isEmpty() ) {
            if ( _phylogeny.getNumberOfExternalNodes() <= LIMIT_FOR_HQ_RENDERING ) {
                _rendering_hints.put( RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY );
            }
            else {
                _rendering_hints.put( RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED );
            }
        }
        if ( getMainPanel().getOptions().isAntialiasScreen() ) {
            if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) {
                _rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF );
            }
            else {
                _rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
            }
            try {
                _rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING,
                                      RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB );
            }
            catch ( final Throwable e ) {
                _rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            }
        }
        else {
            _rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF );
            _rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF );
        }
    }

    final void setTinyFonts() {
        getTreeFontSet().tinyFonts();
    }

    /**
     * Set a phylogeny tree.
     * 
     * @param t
     *            an instance of a Phylogeny
     */
    final void setTree( final Phylogeny t ) {
        _phylogeny = t;
    }

    final void setTreeFile( final File treefile ) {
        _treefile = treefile;
    }

    final private void setUpUrtFactor() {
        final int d = getVisibleRect().width < getVisibleRect().height ? getVisibleRect().width
                : getVisibleRect().height;
        if ( isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram() ) {
            setUrtFactor( ( float ) ( d / ( 2 * getMaxDistanceToRoot() ) ) );
        }
        else {
            final int max_depth = _circ_max_depth;
            if ( max_depth > 0 ) {
                setUrtFactor( d / ( 2 * max_depth ) );
            }
            else {
                setUrtFactor( d / 2 );
            }
        }
        setUrtFactorOv( getUrtFactor() );
    }

    final private void setUrtFactor( final float urt_factor ) {
        _urt_factor = urt_factor;
    }

    final private void setUrtFactorOv( final float urt_factor_ov ) {
        _urt_factor_ov = urt_factor_ov;
    }

    final void setWaitCursor() {
        setCursor( WAIT_CURSOR );
        repaint();
    }

    final void setXcorrectionFactor( final float f ) {
        _x_correction_factor = f;
    }

    final void setXdistance( final float x ) {
        _x_distance = x;
    }

    final void setYdistance( final float y ) {
        _y_distance = y;
    }

    final private void showNodeDataPopup( final MouseEvent e, final PhylogenyNode node ) {
        try {
            if ( ( node.getNodeName().length() > 0 )
                    || ( node.getNodeData().isHasTaxonomy() && !isTaxonomyEmpty( node.getNodeData().getTaxonomy() ) )
                    || ( node.getNodeData().isHasSequence() && !isSequenceEmpty( node.getNodeData().getSequence() ) )
                    || ( node.getNodeData().isHasDate() ) || ( node.getNodeData().isHasDistribution() )
                    || node.getBranchData().isHasConfidences() ) {
                _popup_buffer.setLength( 0 );
                short lines = 0;
                if ( node.getNodeName().length() > 0 ) {
                    lines++;
                    _popup_buffer.append( node.getNodeName() );
                }
                if ( node.getNodeData().isHasTaxonomy() && !isTaxonomyEmpty( node.getNodeData().getTaxonomy() ) ) {
                    lines++;
                    boolean enc_data = false;
                    final Taxonomy tax = node.getNodeData().getTaxonomy();
                    if ( _popup_buffer.length() > 0 ) {
                        _popup_buffer.append( "\n" );
                    }
                    if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                        _popup_buffer.append( "[" );
                        _popup_buffer.append( tax.getTaxonomyCode() );
                        _popup_buffer.append( "]" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " " );
                        }
                        _popup_buffer.append( tax.getScientificName() );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " (" );
                        }
                        else {
                            _popup_buffer.append( "(" );
                        }
                        _popup_buffer.append( tax.getCommonName() );
                        _popup_buffer.append( ")" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( tax.getAuthority() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " (" );
                        }
                        else {
                            _popup_buffer.append( "(" );
                        }
                        _popup_buffer.append( tax.getAuthority() );
                        _popup_buffer.append( ")" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( tax.getRank() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " {" );
                        }
                        else {
                            _popup_buffer.append( "{" );
                        }
                        _popup_buffer.append( tax.getRank() );
                        _popup_buffer.append( "}" );
                        enc_data = true;
                    }
                    if ( tax.getSynonyms().size() > 0 ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " " );
                        }
                        _popup_buffer.append( "[" );
                        int counter = 1;
                        for( final String syn : tax.getSynonyms() ) {
                            if ( !ForesterUtil.isEmpty( syn ) ) {
                                _popup_buffer.append( syn );
                                if ( counter < tax.getSynonyms().size() ) {
                                    _popup_buffer.append( ", " );
                                }
                            }
                            counter++;
                        }
                        _popup_buffer.append( "]" );
                    }
                }
                if ( node.getNodeData().isHasSequence() && !isSequenceEmpty( node.getNodeData().getSequence() ) ) {
                    lines++;
                    boolean enc_data = false;
                    if ( _popup_buffer.length() > 0 ) {
                        _popup_buffer.append( "\n" );
                    }
                    final Sequence seq = node.getNodeData().getSequence();
                    if ( seq.getAccession() != null ) {
                        _popup_buffer.append( "[" );
                        if ( !ForesterUtil.isEmpty( seq.getAccession().getSource() ) ) {
                            _popup_buffer.append( seq.getAccession().getSource() );
                            _popup_buffer.append( "=" );
                        }
                        _popup_buffer.append( seq.getAccession().getValue() );
                        _popup_buffer.append( "]" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( seq.getSymbol() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " [" );
                        }
                        else {
                            _popup_buffer.append( "[" );
                        }
                        _popup_buffer.append( seq.getSymbol() );
                        _popup_buffer.append( "]" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( seq.getName() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " " );
                        }
                        _popup_buffer.append( seq.getName() );
                    }
                }
                if ( node.getNodeData().isHasDate() ) {
                    lines++;
                    if ( _popup_buffer.length() > 0 ) {
                        _popup_buffer.append( "\n" );
                    }
                    _popup_buffer.append( node.getNodeData().getDate().asSimpleText() );
                }
                if ( node.getNodeData().isHasDistribution() ) {
                    lines++;
                    if ( _popup_buffer.length() > 0 ) {
                        _popup_buffer.append( "\n" );
                    }
                    _popup_buffer.append( node.getNodeData().getDistribution().asSimpleText() );
                }
                if ( node.getBranchData().isHasConfidences() ) {
                    final List<Confidence> confs = node.getBranchData().getConfidences();
                    for( final Confidence confidence : confs ) {
                        lines++;
                        if ( _popup_buffer.length() > 0 ) {
                            _popup_buffer.append( "\n" );
                        }
                        if ( !ForesterUtil.isEmpty( confidence.getType() ) ) {
                            _popup_buffer.append( "[" );
                            _popup_buffer.append( confidence.getType() );
                            _popup_buffer.append( "] " );
                        }
                        else {
                            _popup_buffer.append( "[?] " );
                        }
                        _popup_buffer.append( FORMATTER_CONFIDENCE.format( ForesterUtil
                                .round( confidence.getValue(), getOptions()
                                        .getNumberOfDigitsAfterCommaForConfidenceValues() ) ) );
                    }
                }
                if ( _popup_buffer.length() > 0 ) {
                    if ( !getConfiguration().isUseNativeUI() ) {
                        _rollover_popup
                                .setBorder( BorderFactory.createLineBorder( getTreeColorSet().getBranchColor() ) );
                        _rollover_popup.setBackground( getTreeColorSet().getBackgroundColor() );
                        if ( isInFoundNodes( node ) ) {
                            _rollover_popup.setForeground( getTreeColorSet().getFoundColor() );
                        }
                        else if ( getControlPanel().isColorAccordingToTaxonomy() ) {
                            _rollover_popup.setForeground( getTaxonomyBasedColor( node ) );
                        }
                        else {
                            _rollover_popup.setForeground( getTreeColorSet().getSequenceColor() );
                        }
                    }
                    else {
                        _rollover_popup.setBorder( BorderFactory.createLineBorder( Color.BLACK ) );
                    }
                    _rollover_popup.setText( _popup_buffer.toString() );
                    _node_desc_popup = PopupFactory.getSharedInstance().getPopup( null,
                                                                                  _rollover_popup,
                                                                                  e.getLocationOnScreen().x + 10,
                                                                                  e.getLocationOnScreen().y
                                                                                          - ( lines * 20 ) );
                    _node_desc_popup.show();
                }
            }
        }
        catch ( final Exception ex ) {
            // Do nothing.
        }
    }

    final private void showNodeEditFrame( final PhylogenyNode n ) {
        if ( _node_frame_index < TreePanel.MAX_NODE_FRAMES ) {
            // pop up edit box for single node
            _node_frames[ _node_frame_index ] = new NodeFrame( n, _phylogeny, this, _node_frame_index, "" );
            _node_frame_index++;
        }
        else {
            JOptionPane.showMessageDialog( this, "too many node windows are open" );
        }
    }

    final private void showNodeFrame( final PhylogenyNode n ) {
        if ( _node_frame_index < TreePanel.MAX_NODE_FRAMES ) {
            // pop up edit box for single node
            _node_frames[ _node_frame_index ] = new NodeFrame( n, _phylogeny, this, _node_frame_index );
            _node_frame_index++;
        }
        else {
            JOptionPane.showMessageDialog( this, "too many node windows are open" );
        }
    }

    /**
     * Find a color for this species name.
     * 
     * @param species
     * @return the species color
     */
    final Color getTaxonomyBasedColor( final PhylogenyNode node ) {
        if ( node.getNodeData().getTaxonomy() == null ) {
            // return non-colorized color
            return getTreeColorSet().getTaxonomyColor();
        }
        return calculateTaxonomyBasedColor( node.getNodeData().getTaxonomy() );
    }

    final Color calculateTaxonomyBasedColor( final Taxonomy tax ) {
        String species = tax.getTaxonomyCode();
        if ( ForesterUtil.isEmpty( species ) ) {
            species = tax.getScientificName();
            if ( ForesterUtil.isEmpty( species ) ) {
                species = tax.getCommonName();
            }
        }
        if ( ForesterUtil.isEmpty( species ) ) {
            return getTreeColorSet().getTaxonomyColor();
        }
        // Look in species hash
        Color c = getControlPanel().getSpeciesColors().get( species );
        if ( c == null ) {
            c = Util.calculateColorFromString( species );
            getControlPanel().getSpeciesColors().put( species, c );
        }
        return c;
    }

    final void subTree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot get a sub/super tree in unrooted display",
                                           "Attempt to get sub/super tree in unrooted display",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( node.isExternal() ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot get a subtree of a external node",
                                           "Attempt to get subtree of external node",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( node.isRoot() && ( _subtree_index < 1 ) ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot get a subtree of the root node",
                                           "Attempt to get subtree of root node",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( !node.isExternal() && !node.isRoot() && ( _subtree_index <= ( TreePanel.MAX_SUBTREES - 1 ) ) ) {
            _phylogenies[ _subtree_index++ ] = _phylogeny;
            _phylogeny = _phylogeny.subTree( node );
            updateSubSuperTreeButton();
        }
        else if ( node.isRoot() && ( _subtree_index >= 1 ) ) {
            superTree();
        }
        _main_panel.getControlPanel().showWhole();
        repaint();
    }

    final void superTree() {
        _phylogenies[ _subtree_index ] = null;
        _phylogeny = _phylogenies[ --_subtree_index ];
        updateSubSuperTreeButton();
    }

    final void swap( final PhylogenyNode node ) {
        if ( !node.isExternal() ) {
            _phylogeny.swapChildren( node );
        }
        repaint();
    }

    final private void switchDisplaygetPhylogenyGraphicsType() {
        switch ( getPhylogenyGraphicsType() ) {
            case RECTANGULAR:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE );
                break;
            case EURO_STYLE:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.ROUNDED );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.ROUNDED );
                break;
            case ROUNDED:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CURVED );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CURVED );
                break;
            case CURVED:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR );
                break;
            case TRIANGULAR:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CONVEX );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CONVEX );
                break;
            case CONVEX:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                break;
            case UNROOTED:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                break;
            case CIRCULAR:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                break;
            default:
                throw new IllegalStateException( "unkwnown display type: " + getPhylogenyGraphicsType() );
        }
        if ( getControlPanel().getDynamicallyHideData() != null ) {
            if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
                getControlPanel().getDynamicallyHideData().setEnabled( false );
            }
            else {
                getControlPanel().getDynamicallyHideData().setEnabled( true );
            }
        }
        if ( isPhyHasBranchLengths() && ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
            getControlPanel().setDrawPhylogramEnabled( true );
        }
        else {
            getControlPanel().setDrawPhylogramEnabled( false );
        }
        if ( getMainPanel().getMainFrame() == null ) {
            // Must be "E" applet version.
            ( ( ArchaeopteryxE ) ( ( MainPanelApplets ) getMainPanel() ).getApplet() )
                    .setSelectedTypeInTypeMenu( getPhylogenyGraphicsType() );
        }
        else {
            getMainPanel().getMainFrame().setSelectedTypeInTypeMenu( getPhylogenyGraphicsType() );
        }
    }

    final void taxColor() {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        setWaitCursor();
        Util.colorPhylogenyAccordingToExternalTaxonomy( _phylogeny, this );
        _control_panel.setColorBranches( true );
        if ( _control_panel.getColorBranchesCb() != null ) {
            _control_panel.getColorBranchesCb().setSelected( true );
        }
        setArrowCursor();
        repaint();
    }

    final void updateOvSettings() {
        switch ( getOptions().getOvPlacement() ) {
            case LOWER_LEFT:
                setOvXPosition( OV_BORDER );
                setOvYPosition( ForesterUtil.roundToInt( getVisibleRect().height - OV_BORDER - getOvMaxHeight() ) );
                setOvYStart( ForesterUtil.roundToInt( getOvYPosition() + ( getOvMaxHeight() / 2 ) ) );
                break;
            case LOWER_RIGHT:
                setOvXPosition( ForesterUtil.roundToInt( getVisibleRect().width - OV_BORDER - getOvMaxWidth() ) );
                setOvYPosition( ForesterUtil.roundToInt( getVisibleRect().height - OV_BORDER - getOvMaxHeight() ) );
                setOvYStart( ForesterUtil.roundToInt( getOvYPosition() + ( getOvMaxHeight() / 2 ) ) );
                break;
            case UPPER_RIGHT:
                setOvXPosition( ForesterUtil.roundToInt( getVisibleRect().width - OV_BORDER - getOvMaxWidth() ) );
                setOvYPosition( OV_BORDER );
                setOvYStart( ForesterUtil.roundToInt( OV_BORDER + ( getOvMaxHeight() / 2 ) ) );
                break;
            default:
                setOvXPosition( OV_BORDER );
                setOvYPosition( OV_BORDER );
                setOvYStart( ForesterUtil.roundToInt( OV_BORDER + ( getOvMaxHeight() / 2 ) ) );
                break;
        }
    }

    final void updateOvSizes() {
        if ( ( getWidth() > 1.05 * getVisibleRect().width ) || ( getHeight() > 1.05 * getVisibleRect().height ) ) {
            setOvOn( true );
            float l = getLongestExtNodeInfo();
            final float w_ratio = getOvMaxWidth() / getWidth();
            l *= w_ratio;
            final int ext_nodes = _phylogeny.getRoot().getNumberOfExternalNodes();
            setOvYDistance( getOvMaxHeight() / ( 2 * ext_nodes ) );
            float ov_xdist = 0;
            if ( !isNonLinedUpCladogram() && !isUniformBranchLengthsForCladogram() ) {
                ov_xdist = ( ( getOvMaxWidth() - l ) / ( ext_nodes ) );
            }
            else {
                ov_xdist = ( ( getOvMaxWidth() - l ) / ( PhylogenyMethods.calculateMaxDepth( _phylogeny ) ) );
            }
            float ydist = ( float ) ( ( getOvMaxWidth() / ( ext_nodes * 2.0 ) ) );
            if ( ov_xdist < 0.0 ) {
                ov_xdist = 0.0f;
            }
            if ( ydist < 0.0 ) {
                ydist = 0.0f;
            }
            setOvXDistance( ov_xdist );
            final double height = _phylogeny.getHeight();
            if ( height > 0 ) {
                final float ov_corr = ( float ) ( ( ( getOvMaxWidth() - l ) - getOvXDistance() ) / height );
                setOvXcorrectionFactor( ov_corr > 0 ? ov_corr : 0 );
            }
            else {
                setOvXcorrectionFactor( 0 );
            }
        }
        else {
            setOvOn( false );
        }
    }

    final void updateSubSuperTreeButton() {
        if ( _subtree_index < 1 ) {
            getControlPanel().deactivateButtonToReturnToSuperTree();
        }
        else {
            getControlPanel().activateButtonToReturnToSuperTree( _subtree_index );
        }
    }

    final void zoomInDomainStructure() {
        if ( _domain_structure_width < 2000 ) {
            _domain_structure_width *= 1.2;
        }
    }

    final void zoomOutDomainStructure() {
        if ( _domain_structure_width > 20 ) {
            _domain_structure_width *= 0.8;
        }
    }

    final private static void drawString( final int i, final double x, final double y, final Graphics2D g ) {
        g.drawString( String.valueOf( i ), ( int ) ( x + 0.5 ), ( int ) ( y + 0.5 ) );
    }

    final private static void drawString( final String str, final double x, final double y, final Graphics2D g ) {
        g.drawString( str, ( int ) ( x + 0.5 ), ( int ) ( y + 0.5 ) );
    }

    final private static boolean isSequenceEmpty( final Sequence seq ) {
        return ( seq.getAccession() == null ) && ForesterUtil.isEmpty( seq.getName() )
                && ForesterUtil.isEmpty( seq.getSymbol() );
    }

    final private static boolean isTaxonomyEmpty( final Taxonomy tax ) {
        return ( ( tax.getIdentifier() == null ) && ForesterUtil.isEmpty( tax.getTaxonomyCode() )
                && ForesterUtil.isEmpty( tax.getCommonName() ) && ForesterUtil.isEmpty( tax.getScientificName() ) && ( tax
                .getSynonyms().isEmpty() ) );
    }

    final private static boolean plusPressed( final int key_code ) {
        return ( ( key_code == KeyEvent.VK_ADD ) || ( key_code == KeyEvent.VK_PLUS )
                || ( key_code == KeyEvent.VK_EQUALS ) || ( key_code == KeyEvent.VK_SEMICOLON ) || ( key_code == KeyEvent.VK_1 ) );
    }

    final private class SubtreeColorizationActionListener implements ActionListener {

        JColorChooser _chooser;
        PhylogenyNode _node;

        SubtreeColorizationActionListener( final JColorChooser chooser, final PhylogenyNode node ) {
            _chooser = chooser;
            _node = node;
        }

        @Override
        public void actionPerformed( final ActionEvent e ) {
            final Color c = _chooser.getColor();
            if ( c != null ) {
                colorizeSubtree( c, _node );
            }
        }
    }
}
