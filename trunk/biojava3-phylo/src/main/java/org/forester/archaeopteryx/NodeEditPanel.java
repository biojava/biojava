// $Id: NodeEditPanel.java,v 1.19 2009/10/30 03:00:51 cmzmasek Exp $
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
// WWW: www.phylosoft.org/

package org.forester.archaeopteryx;

import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.math.BigDecimal;
import java.net.URI;
import java.text.ParseException;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JEditorPane;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTree;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.text.Position;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Reference;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;
import org.forester.util.ForesterUtil;

class NodeEditPanel extends JPanel {

    private static final long                            serialVersionUID = 5120159904388100771L;
    private final JTree                                  _tree;
    private final JEditorPane                            _pane;
    private final PhylogenyNode                          _my_node;
    private final TreePanel                              _tree_panel;
    private final Map<DefaultMutableTreeNode, TagNumber> _map;

    public NodeEditPanel( final PhylogenyNode phylogeny_node, final TreePanel tree_panel ) {
        _map = new HashMap<DefaultMutableTreeNode, TagNumber>();
        _my_node = phylogeny_node;
        _tree_panel = tree_panel;
        String node_name = "";
        if ( !ForesterUtil.isEmpty( phylogeny_node.getNodeName() ) ) {
            node_name = phylogeny_node.getNodeName() + " ";
        }
        final DefaultMutableTreeNode top = new DefaultMutableTreeNode( "Node " + node_name );
        createNodes( top, phylogeny_node );
        _tree = new JTree( top );
        getJTree().setEditable( true );
        getJTree().setFocusable( true );
        getJTree().setToggleClickCount( 1 );
        getJTree().setInvokesStopCellEditing( true );
        final JScrollPane tree_view = new JScrollPane( getJTree() );
        _pane = new JEditorPane();
        _pane.setEditable( true );
        final JScrollPane data_view = new JScrollPane( _pane );
        final JSplitPane split_pane = new JSplitPane( JSplitPane.VERTICAL_SPLIT );
        split_pane.setTopComponent( tree_view );
        // split_pane.setBottomComponent( data_view );
        data_view.setMinimumSize( Constants.NODE_PANEL_SPLIT_MINIMUM_SIZE );
        tree_view.setMinimumSize( Constants.NODE_PANEL_SPLIT_MINIMUM_SIZE );
        // split_pane.setDividerLocation( 400 );
        split_pane.setPreferredSize( Constants.NODE_PANEL_SIZE );
        add( split_pane );
        getJTree().getSelectionModel().setSelectionMode( TreeSelectionModel.SINGLE_TREE_SELECTION );
        getJTree().addKeyListener( new KeyListener() {

            @Override
            public void keyPressed( final KeyEvent e ) {
                keyEvent( e );
            }

            @Override
            public void keyReleased( final KeyEvent e ) {
                keyEvent( e );
            }

            @Override
            public void keyTyped( final KeyEvent e ) {
                keyEvent( e );
            }
        } );
        for( int i = 0; i < getJTree().getRowCount(); i++ ) {
            getJTree().expandRow( i );
        }
        collapsePath( NodePanel.BASIC );
        collapsePath( NodePanel.TAXONOMY );
        collapsePath( NodePanel.SEQUENCE );
        collapsePath( NodePanel.EVENTS );
        collapsePath( NodePanel.DATE );
        collapsePath( NodePanel.DISTRIBUTION );
        collapsePath( NodePanel.LIT_REFERENCE );
        getJTree().addTreeSelectionListener( new TreeSelectionListener() {

            @Override
            public void valueChanged( final TreeSelectionEvent e ) {
                final TreePath new_path = e.getNewLeadSelectionPath();
                final TreePath old_path = e.getOldLeadSelectionPath();
                if ( new_path != null ) {
                    writeBack( ( DefaultMutableTreeNode ) new_path.getLastPathComponent() );
                }
                if ( old_path != null ) {
                    writeBack( ( DefaultMutableTreeNode ) old_path.getLastPathComponent() );
                }
            }
        } );
    }

    private void addBasics( final DefaultMutableTreeNode top, final PhylogenyNode phylogeny_node, final String name ) {
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelementEditable( category, NodePanel.NODE_NAME, phylogeny_node.getNodeName(), PHYLOXML_TAG.NODE_NAME );
        String bl = "";
        if ( phylogeny_node.getDistanceToParent() != PhylogenyNode.DISTANCE_DEFAULT ) {
            bl = ForesterUtil.FORMATTER_6.format( phylogeny_node.getDistanceToParent() );
        }
        addSubelementEditable( category, NodePanel.NODE_BRANCH_LENGTH, bl, PHYLOXML_TAG.NODE_BRANCH_LENGTH );
        int counter = 0;
        if ( phylogeny_node.getBranchData().isHasConfidences() ) {
            for( int i = phylogeny_node.getBranchData().getConfidences().size() - 1; i >= 0; i-- ) {
                if ( phylogeny_node.getBranchData().getConfidences().get( i ).getValue() == Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                    phylogeny_node.getBranchData().getConfidences().remove( i );
                }
            }
            for( final PhylogenyData conf : phylogeny_node.getBranchData().getConfidences() ) {
                final Confidence my_conf = ( Confidence ) ( conf );
                addSubelementEditable( category,
                                       NodePanel.CONFIDENCE + " [" + counter + "]",
                                       ForesterUtil.FORMATTER_6.format( my_conf.getValue() ),
                                       PHYLOXML_TAG.CONFIDENCE_VALUE,
                                       NodePanel.CONFIDENCE_TYPE,
                                       my_conf.getType(),
                                       PHYLOXML_TAG.CONFIDENCE_TYPE,
                                       counter++ );
            }
        }
        addSubelementEditable( category,
                               NodePanel.CONFIDENCE + " [" + counter + "]",
                               "",
                               PHYLOXML_TAG.CONFIDENCE_VALUE,
                               NodePanel.CONFIDENCE_TYPE,
                               "",
                               PHYLOXML_TAG.CONFIDENCE_TYPE,
                               counter );
    }

    //    private void addAnnotation( final DefaultMutableTreeNode top, final Annotation ann, final String name ) {
    //        DefaultMutableTreeNode category;
    //        category = new DefaultMutableTreeNode( name );
    //        top.add( category );
    //        addSubelementEditable( category, "Reference", ann.getRef() , PHYLOXML_TAG.);
    //        addSubelementEditable( category, "Description", ann.getDesc() , PHYLOXML_TAG.);
    //        addSubelementEditable( category, "Source", ann.getSource(), PHYLOXML_TAG. );
    //        addSubelementEditable( category, "Type", ann.getType(), PHYLOXML_TAG. );
    //        addSubelementEditable( category, "Evidence", ann.getEvidence() , PHYLOXML_TAG.);
    //        if ( ann.getConfidence() != null ) {
    //            addSubelementEditable( category, "Confidence", ann.getConfidence().asText().toString() , PHYLOXML_TAG.);
    //        }
    //        if ( ann.getProperties() != null ) {
    //            addProperties( category, ann.getProperties(), "Properties", PHYLOXML_TAG. );
    //        }
    //    }
    //    private void addAnnotations( final DefaultMutableTreeNode top,
    //                                 final List<PhylogenyData> annotations,
    //                                 final DefaultMutableTreeNode category ) {
    //        if ( ( annotations != null ) && ( annotations.size() > 0 ) ) {
    //            category.add( new DefaultMutableTreeNode( "Annotations" ) );
    //            final DefaultMutableTreeNode last = top.getLastLeaf();
    //            int i = 0;
    //            for( final PhylogenyData ann : annotations ) {
    //                addAnnotation( last, ( Annotation ) ann, "Annotation " + ( i++ ) );
    //            }
    //        }
    //    }
    private void addDate( final DefaultMutableTreeNode top, Date date, final String name ) {
        if ( date == null ) {
            date = new Date();
        }
        DefaultMutableTreeNode category;
        category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelementEditable( category, NodePanel.DATE_DESCRIPTION, date.getDesc(), PHYLOXML_TAG.DATE_DESCRIPTION );
        addSubelementEditable( category, NodePanel.DATE_VALUE, String.valueOf( date.getValue() != null ? date
                .getValue() : "" ), PHYLOXML_TAG.DATE_VALUE );
        addSubelementEditable( category, NodePanel.DATE_MIN, String
                .valueOf( date.getMin() != null ? date.getMin() : "" ), PHYLOXML_TAG.DATE_MIN );
        addSubelementEditable( category, NodePanel.DATE_MAX, String
                .valueOf( date.getMax() != null ? date.getMax() : "" ), PHYLOXML_TAG.DATE_MAX );
        addSubelementEditable( category, NodePanel.DATE_UNIT, date.getUnit(), PHYLOXML_TAG.DATE_UNIT );
    }

    private void addDistribution( final DefaultMutableTreeNode top, Distribution dist, final String name ) {
        if ( dist == null ) {
            dist = new Distribution( "" );
        }
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelementEditable( category, NodePanel.DIST_DESCRIPTION, dist.getDesc(), PHYLOXML_TAG.DIST_DESC );
        addSubelementEditable( category,
                               NodePanel.DIST_GEODETIC_DATUM,
                               dist.getGeodeticDatum(),
                               PHYLOXML_TAG.DIST_GEODETIC );
        addSubelementEditable( category, NodePanel.DIST_LATITUDE, String.valueOf( dist.getLatitude() != null ? dist
                .getLatitude() : "" ), PHYLOXML_TAG.DIST_LAT );
        addSubelementEditable( category, NodePanel.DIST_LONGITUDE, String.valueOf( dist.getLongitude() != null ? dist
                .getLongitude() : "" ), PHYLOXML_TAG.DIST_LONG );
        addSubelementEditable( category, NodePanel.DIST_ALTITUDE, String.valueOf( dist.getAltitude() != null ? dist
                .getAltitude() : "" ), PHYLOXML_TAG.DIST_ALT );
    }

    private void addEvents( final DefaultMutableTreeNode top, Event events, final String name ) {
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        if ( events == null ) {
            events = new Event();
        }
        top.add( category );
        addSubelementEditable( category,
                               NodePanel.EVENTS_DUPLICATIONS,
                               String.valueOf( events.getNumberOfDuplications() >= 0 ? events.getNumberOfDuplications()
                                       : 0 ),
                               PHYLOXML_TAG.EVENTS_DUPLICATIONS );
        addSubelementEditable( category,
                               NodePanel.EVENTS_SPECIATIONS,
                               String.valueOf( events.getNumberOfSpeciations() >= 0 ? events.getNumberOfSpeciations()
                                       : 0 ),
                               PHYLOXML_TAG.EVENTS_SPECIATIONS );
        addSubelementEditable( category,
                               NodePanel.EVENTS_GENE_LOSSES,
                               String
                                       .valueOf( events.getNumberOfGeneLosses() >= 0 ? events.getNumberOfGeneLosses()
                                               : 0 ),
                               PHYLOXML_TAG.EVENTS_GENE_LOSSES );
    }

    private void addMapping( final DefaultMutableTreeNode mtn, final TagNumber tag ) {
        if ( getMap().containsKey( mtn ) ) {
            throw new IllegalArgumentException( "key " + mtn + " already present" );
        }
        if ( getMap().containsValue( tag ) ) {
            throw new IllegalArgumentException( "value " + tag + " already present" );
        }
        getMap().put( mtn, tag );
    }

    private void addReference( final DefaultMutableTreeNode top, Reference ref, final String name ) {
        if ( ref == null ) {
            ref = new Reference( "" );
        }
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        addSubelementEditable( category, NodePanel.LIT_REFERENCE_DESC, ref.getValue(), PHYLOXML_TAG.LIT_REFERENCE_DESC );
        addSubelementEditable( category, NodePanel.LIT_REFERENCE_DOI, ref.getDoi(), PHYLOXML_TAG.LIT_REFERENCE_DOI );
    }

    private void addSequence( final DefaultMutableTreeNode top, Sequence seq, final String name ) {
        if ( seq == null ) {
            seq = new Sequence();
        }
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        Accession acc = seq.getAccession();
        if ( acc == null ) {
            acc = new Accession( "", "" );
        }
        addSubelementEditable( category, NodePanel.SEQ_NAME, seq.getName(), PHYLOXML_TAG.SEQ_NAME );
        addSubelementEditable( category, NodePanel.SEQ_SYMBOL, seq.getSymbol(), PHYLOXML_TAG.SEQ_SYMBOL );
        addSubelementEditable( category,
                               NodePanel.SEQ_ACCESSION,
                               acc.getValue(),
                               PHYLOXML_TAG.SEQ_ACC_VALUE,
                               "Source",
                               acc.getSource(),
                               PHYLOXML_TAG.SEQ_ACC_SOURCE );
        addSubelementEditable( category, NodePanel.SEQ_LOCATION, seq.getLocation(), PHYLOXML_TAG.SEQ_LOCATION );
        addSubelementEditable( category, NodePanel.SEQ_TYPE, seq.getType(), PHYLOXML_TAG.SEQ_TYPE );
        addSubelementEditable( category, NodePanel.SEQ_MOL_SEQ, seq.getMolecularSequence(), PHYLOXML_TAG.SEQ_MOL_SEQ );
        if ( seq.getUri() != null ) {
            addSubelementEditable( category, NodePanel.SEQ_URI, seq.getUri().toString(), PHYLOXML_TAG.SEQ_URI );
        }
        else {
            addSubelementEditable( category, NodePanel.SEQ_URI, "", PHYLOXML_TAG.SEQ_URI );
        }
        //  addAnnotations( top, seq.getAnnotations(), category );
    }

    private void addSubelementEditable( final DefaultMutableTreeNode node,
                                        final String name,
                                        final String value,
                                        final PHYLOXML_TAG phyloxml_tag ) {
        addSubelementEditable( node, name, value, phyloxml_tag, 0 );
    }

    private void addSubelementEditable( final DefaultMutableTreeNode node,
                                        final String name,
                                        final String value,
                                        final PHYLOXML_TAG phyloxml_tag,
                                        final int number ) {
        String my_value = value;
        if ( ForesterUtil.isEmpty( my_value ) ) {
            my_value = "";
        }
        final DefaultMutableTreeNode name_node = new DefaultMutableTreeNode( name );
        final DefaultMutableTreeNode value_node = new DefaultMutableTreeNode( my_value );
        name_node.add( value_node );
        node.add( name_node );
        addMapping( name_node, new TagNumber( phyloxml_tag, number ) );
    }

    private void addSubelementEditable( final DefaultMutableTreeNode node,
                                        final String name,
                                        final String value,
                                        final PHYLOXML_TAG phyloxml_value_tag,
                                        final String source_name,
                                        final String source_value,
                                        final PHYLOXML_TAG phyloxml_source_tag ) {
        addSubelementEditable( node, name, value, phyloxml_value_tag, source_name, source_value, phyloxml_source_tag, 0 );
    }

    private void addSubelementEditable( final DefaultMutableTreeNode node,
                                        final String name,
                                        final String value,
                                        final PHYLOXML_TAG phyloxml_value_tag,
                                        final String source_name,
                                        final String source_value,
                                        final PHYLOXML_TAG phyloxml_source_tag,
                                        final int number ) {
        String my_value = value;
        if ( ForesterUtil.isEmpty( my_value ) ) {
            my_value = "";
        }
        String my_source_value = source_value;
        if ( ForesterUtil.isEmpty( my_source_value ) ) {
            my_source_value = "";
        }
        final DefaultMutableTreeNode name_node = new DefaultMutableTreeNode( name );
        final DefaultMutableTreeNode source_name_node = new DefaultMutableTreeNode( source_name );
        final DefaultMutableTreeNode source_value_node = new DefaultMutableTreeNode( my_source_value );
        final DefaultMutableTreeNode value_node = new DefaultMutableTreeNode( my_value );
        name_node.add( source_name_node );
        source_name_node.add( source_value_node );
        name_node.add( value_node );
        node.add( name_node );
        addMapping( name_node, new TagNumber( phyloxml_value_tag, number ) );
        addMapping( source_name_node, new TagNumber( phyloxml_source_tag, number ) );
    }

    private void addTaxonomy( final DefaultMutableTreeNode top, Taxonomy tax, final String name ) {
        if ( tax == null ) {
            tax = new Taxonomy();
        }
        final DefaultMutableTreeNode category = new DefaultMutableTreeNode( name );
        top.add( category );
        Identifier id = tax.getIdentifier();
        if ( id == null ) {
            id = new Identifier();
        }
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_IDENTIFIER,
                               id.getValue(),
                               PHYLOXML_TAG.TAXONOMY_ID_VALUE,
                               "Provider",
                               id.getProvider(),
                               PHYLOXML_TAG.TAXONOMY_ID_PROVIDER );
        addSubelementEditable( category, NodePanel.TAXONOMY_CODE, tax.getTaxonomyCode(), PHYLOXML_TAG.TAXONOMY_CODE );
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_SCIENTIFIC_NAME,
                               tax.getScientificName(),
                               PHYLOXML_TAG.TAXONOMY_SCIENTIFIC_NAME );
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_AUTHORITY,
                               tax.getAuthority(),
                               PHYLOXML_TAG.TAXONOMY_AUTHORITY );
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_COMMON_NAME,
                               tax.getCommonName(),
                               PHYLOXML_TAG.TAXONOMY_COMMON_NAME );
        for( int i = tax.getSynonyms().size() - 1; i >= 0; i-- ) {
            if ( ForesterUtil.isEmpty( tax.getSynonyms().get( i ) ) ) {
                tax.getSynonyms().remove( i );
            }
        }
        int counter = 0;
        for( final String syn : tax.getSynonyms() ) {
            addSubelementEditable( category,
                                   NodePanel.TAXONOMY_SYNONYM + " [" + counter + "]",
                                   syn,
                                   PHYLOXML_TAG.TAXONOMY_SYNONYM,
                                   counter++ );
        }
        addSubelementEditable( category,
                               NodePanel.TAXONOMY_SYNONYM + " [" + counter + "]",
                               "",
                               PHYLOXML_TAG.TAXONOMY_SYNONYM,
                               counter );
        addSubelementEditable( category, NodePanel.TAXONOMY_RANK, tax.getRank(), PHYLOXML_TAG.TAXONOMY_RANK );
        if ( tax.getUri() != null ) {
            addSubelementEditable( category, NodePanel.TAXONOMY_URI, tax.getUri().toString(), PHYLOXML_TAG.TAXONOMY_URI );
        }
        else {
            addSubelementEditable( category, NodePanel.TAXONOMY_URI, "", PHYLOXML_TAG.TAXONOMY_URI );
        }
    }

    private void collapsePath( final String name ) {
        final TreePath tp = getJTree().getNextMatch( name, 0, Position.Bias.Forward );
        if ( tp != null ) {
            getJTree().collapsePath( tp );
        }
    }

    private void createNodes( final DefaultMutableTreeNode top, final PhylogenyNode phylogeny_node ) {
        addBasics( top, phylogeny_node, NodePanel.BASIC );
        addTaxonomy( top, phylogeny_node.getNodeData().getTaxonomy(), NodePanel.TAXONOMY );
        addSequence( top, phylogeny_node.getNodeData().getSequence(), NodePanel.SEQUENCE );
        if ( !phylogeny_node.isExternal() ) {
            addEvents( top, phylogeny_node.getNodeData().getEvent(), NodePanel.EVENTS );
        }
        addDate( top, phylogeny_node.getNodeData().getDate(), NodePanel.DATE );
        addDistribution( top, phylogeny_node.getNodeData().getDistribution(), NodePanel.DISTRIBUTION );
        addReference( top, phylogeny_node.getNodeData().getReference(), NodePanel.LIT_REFERENCE );
        //  addProperties( top, phylogeny_node.getNodeData().getProperties(), "Properties" );
    }

    private void formatError( final DefaultMutableTreeNode mtn, final PhyloXmlDataFormatException e ) {
        JOptionPane.showMessageDialog( this, e.getMessage(), "Format error", JOptionPane.ERROR_MESSAGE );
        mtn.setUserObject( "" );
        getJTree().repaint();
    }

    private JTree getJTree() {
        return _tree;
    }

    private Map<DefaultMutableTreeNode, TagNumber> getMap() {
        return _map;
    }

    private TagNumber getMapping( final DefaultMutableTreeNode mtn ) {
        return getMap().get( mtn );
    }

    PhylogenyNode getMyNode() {
        return _my_node;
    }

    private DefaultMutableTreeNode getSelectedTreeNode() {
        final TreePath selectionPath = getJTree().getSelectionPath();
        if ( selectionPath != null ) {
            final Object[] path = selectionPath.getPath();
            if ( path.length > 0 ) {
                return ( DefaultMutableTreeNode ) path[ path.length - 1 ]; // Last node
            }
        }
        return null;
    }

    private TreePanel getTreePanel() {
        return _tree_panel;
    }

    private void keyEvent( final KeyEvent e ) {
        if ( e.getKeyCode() == KeyEvent.VK_ENTER ) {
            writeBack( getSelectedTreeNode() );
        }
    }

    private BigDecimal parseBigDecimal( final DefaultMutableTreeNode mtn, final String value ) {
        if ( ForesterUtil.isEmpty( value ) ) {
            return new BigDecimal( 0 );
        }
        BigDecimal i = null;
        try {
            i = new BigDecimal( value );
        }
        catch ( final NumberFormatException e ) {
            JOptionPane.showMessageDialog( this, "illegal value: " + value, "Error", JOptionPane.ERROR_MESSAGE );
            mtn.setUserObject( "" );
        }
        return i;
    }

    private int parsePositiveInt( final DefaultMutableTreeNode mtn, final String value ) {
        if ( ForesterUtil.isEmpty( value ) ) {
            return 0;
        }
        int i = -1;
        try {
            i = ForesterUtil.parseInt( value );
        }
        catch ( final ParseException e ) {
            JOptionPane.showMessageDialog( this, "illegal value: " + value, "Error", JOptionPane.ERROR_MESSAGE );
            mtn.setUserObject( "" );
        }
        if ( i < 0 ) {
            JOptionPane.showMessageDialog( this, "illegal value: " + value, "Error", JOptionPane.ERROR_MESSAGE );
            mtn.setUserObject( "" );
        }
        return i;
    }

    void writeAll() {
        for( int i = 0; i < getJTree().getRowCount(); i++ ) {
            final TreePath p = getJTree().getPathForRow( i );
            writeBack( ( DefaultMutableTreeNode ) p.getLastPathComponent() );
        }
    }

    private void writeBack( final DefaultMutableTreeNode mtn ) {
        if ( !getMap().containsKey( mtn ) ) {
            final DefaultMutableTreeNode parent = ( DefaultMutableTreeNode ) mtn.getParent();
            if ( getMap().containsKey( parent ) ) {
                writeBack( mtn, getMapping( parent ) );
            }
        }
    }

    private void writeBack( final DefaultMutableTreeNode mtn, final TagNumber tag_number ) {
        if ( tag_number == null ) {
            return;
        }
        String value = mtn.toString();
        if ( value == null ) {
            value = "";
        }
        value = value.replaceAll( "\\s+", " " );
        value = value.trim();
        mtn.setUserObject( value );
        getJTree().repaint();
        final PHYLOXML_TAG tag = tag_number.getTag();
        final int number = tag_number.getNumber();
        switch ( tag ) {
            case NODE_NAME:
                getMyNode().setName( value );
                break;
            case NODE_BRANCH_LENGTH:
                if ( ForesterUtil.isEmpty( value ) ) {
                    getMyNode().setDistanceToParent( PhylogenyNode.DISTANCE_DEFAULT );
                }
                else {
                    try {
                        getMyNode().setDistanceToParent( ForesterUtil.parseDouble( value ) );
                    }
                    catch ( final ParseException e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "failed to parse branch length from: " + value,
                                                       "Error",
                                                       JOptionPane.ERROR_MESSAGE );
                        mtn.setUserObject( "" );
                    }
                }
                break;
            case CONFIDENCE_VALUE:
                double confidence = Confidence.CONFIDENCE_DEFAULT_VALUE;
                if ( !ForesterUtil.isEmpty( value ) ) {
                    try {
                        confidence = ForesterUtil.parseDouble( value );
                    }
                    catch ( final ParseException e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "failed to parse confidence value from: " + value,
                                                       "Error",
                                                       JOptionPane.ERROR_MESSAGE );
                        mtn.setUserObject( "" );
                        break;
                    }
                }
                if ( getMyNode().getBranchData().getConfidences().size() < number ) {
                    throw new IllegalStateException();
                }
                else if ( getMyNode().getBranchData().getConfidences().size() == number ) {
                    if ( confidence >= 0 ) {
                        getMyNode().getBranchData().getConfidences().add( new Confidence( confidence, "unknown" ) );
                    }
                }
                else {
                    final String type = getMyNode().getBranchData().getConfidences().get( number ).getType();
                    getMyNode().getBranchData().getConfidences().set( number, new Confidence( confidence, type ) );
                }
                break;
            case CONFIDENCE_TYPE:
                if ( getMyNode().getBranchData().getConfidences().size() < number ) {
                    throw new IllegalStateException();
                }
                else if ( getMyNode().getBranchData().getConfidences().size() == number ) {
                    if ( !ForesterUtil.isEmpty( value ) ) {
                        getMyNode().getBranchData().getConfidences().add( new Confidence( 0, value ) );
                    }
                }
                else {
                    final double v = getMyNode().getBranchData().getConfidences().get( number ).getValue();
                    getMyNode().getBranchData().getConfidences().set( number, new Confidence( v, value ) );
                }
                break;
            case TAXONOMY_CODE:
                ForesterUtil.ensurePresenceOfTaxonomy( getMyNode() );
                try {
                    getMyNode().getNodeData().getTaxonomy().setTaxonomyCode( value );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case TAXONOMY_SCIENTIFIC_NAME:
                ForesterUtil.ensurePresenceOfTaxonomy( getMyNode() );
                getMyNode().getNodeData().getTaxonomy().setScientificName( value );
                break;
            case TAXONOMY_COMMON_NAME:
                ForesterUtil.ensurePresenceOfTaxonomy( getMyNode() );
                getMyNode().getNodeData().getTaxonomy().setCommonName( value );
                break;
            case TAXONOMY_RANK:
                ForesterUtil.ensurePresenceOfTaxonomy( getMyNode() );
                try {
                    getMyNode().getNodeData().getTaxonomy().setRank( value.toLowerCase() );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case TAXONOMY_AUTHORITY:
                ForesterUtil.ensurePresenceOfTaxonomy( getMyNode() );
                getMyNode().getNodeData().getTaxonomy().setAuthority( value );
                break;
            case TAXONOMY_URI:
                ForesterUtil.ensurePresenceOfTaxonomy( getMyNode() );
                if ( ForesterUtil.isEmpty( value ) ) {
                    getMyNode().getNodeData().getTaxonomy().setUri( null );
                }
                else {
                    try {
                        getMyNode().getNodeData().getTaxonomy().setUri( new Uri( URI.create( value ) ) );
                    }
                    catch ( final Exception e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "failed to parse URI from: " + value,
                                                       "Error",
                                                       JOptionPane.ERROR_MESSAGE );
                        mtn.setUserObject( "" );
                    }
                }
                break;
            case TAXONOMY_SYNONYM:
                if ( getMyNode().getNodeData().getTaxonomy().getSynonyms().size() < number ) {
                    throw new IllegalStateException();
                }
                else if ( getMyNode().getNodeData().getTaxonomy().getSynonyms().size() == number ) {
                    if ( !ForesterUtil.isEmpty( value ) ) {
                        ForesterUtil.ensurePresenceOfTaxonomy( getMyNode() );
                        getMyNode().getNodeData().getTaxonomy().getSynonyms().add( value );
                    }
                }
                else {
                    getMyNode().getNodeData().getTaxonomy().getSynonyms().set( number, value );
                }
                break;
            case TAXONOMY_ID_VALUE:
                ForesterUtil.ensurePresenceOfTaxonomy( getMyNode() );
                if ( getMyNode().getNodeData().getTaxonomy().getIdentifier() == null ) {
                    getMyNode().getNodeData().getTaxonomy().setIdentifier( new Identifier( value ) );
                }
                else {
                    final String provider = getMyNode().getNodeData().getTaxonomy().getIdentifier().getProvider();
                    getMyNode().getNodeData().getTaxonomy().setIdentifier( new Identifier( value, provider ) );
                }
                break;
            case TAXONOMY_ID_PROVIDER:
                ForesterUtil.ensurePresenceOfTaxonomy( getMyNode() );
                if ( getMyNode().getNodeData().getTaxonomy().getIdentifier() == null ) {
                    getMyNode().getNodeData().getTaxonomy().setIdentifier( new Identifier( "", value ) );
                }
                else {
                    final String v = getMyNode().getNodeData().getTaxonomy().getIdentifier().getValue();
                    getMyNode().getNodeData().getTaxonomy().setIdentifier( new Identifier( v, value ) );
                }
                break;
            case SEQ_LOCATION:
                ForesterUtil.ensurePresenceOfSequence( getMyNode() );
                getMyNode().getNodeData().getSequence().setLocation( value );
                break;
            case SEQ_MOL_SEQ:
                ForesterUtil.ensurePresenceOfSequence( getMyNode() );
                getMyNode().getNodeData().getSequence().setMolecularSequence( value );
                break;
            case SEQ_NAME:
                ForesterUtil.ensurePresenceOfSequence( getMyNode() );
                getMyNode().getNodeData().getSequence().setName( value );
                break;
            case SEQ_SYMBOL:
                ForesterUtil.ensurePresenceOfSequence( getMyNode() );
                try {
                    getMyNode().getNodeData().getSequence().setSymbol( value );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case SEQ_TYPE:
                ForesterUtil.ensurePresenceOfSequence( getMyNode() );
                try {
                    getMyNode().getNodeData().getSequence().setType( value.toLowerCase() );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case SEQ_ACC_SOURCE:
                ForesterUtil.ensurePresenceOfSequence( getMyNode() );
                if ( getMyNode().getNodeData().getSequence().getAccession() == null ) {
                    getMyNode().getNodeData().getSequence().setAccession( new Accession( "", value ) );
                }
                else {
                    final String v = getMyNode().getNodeData().getSequence().getAccession().getValue();
                    getMyNode().getNodeData().getSequence().setAccession( new Accession( v, value ) );
                }
                break;
            case SEQ_ACC_VALUE:
                ForesterUtil.ensurePresenceOfSequence( getMyNode() );
                if ( getMyNode().getNodeData().getSequence().getAccession() == null ) {
                    getMyNode().getNodeData().getSequence().setAccession( new Accession( value, "" ) );
                }
                else {
                    final String source = getMyNode().getNodeData().getSequence().getAccession().getSource();
                    getMyNode().getNodeData().getSequence().setAccession( new Accession( value, source ) );
                }
                break;
            case SEQ_URI:
                ForesterUtil.ensurePresenceOfSequence( getMyNode() );
                if ( ForesterUtil.isEmpty( value ) ) {
                    getMyNode().getNodeData().getSequence().setUri( null );
                }
                else {
                    try {
                        getMyNode().getNodeData().getSequence().setUri( new Uri( URI.create( value ) ) );
                    }
                    catch ( final Exception e ) {
                        JOptionPane.showMessageDialog( this,
                                                       "failed to parse URI from: " + value,
                                                       "Error",
                                                       JOptionPane.ERROR_MESSAGE );
                        mtn.setUserObject( "" );
                    }
                }
                break;
            case LIT_REFERENCE_DESC:
                if ( !getMyNode().getNodeData().isHasReference() ) {
                    getMyNode().getNodeData().setReference( new Reference( "" ) );
                }
                getMyNode().getNodeData().getReference().setValue( value );
                break;
            case LIT_REFERENCE_DOI:
                if ( !getMyNode().getNodeData().isHasReference() ) {
                    getMyNode().getNodeData().setReference( new Reference( "" ) );
                }
                try {
                    getMyNode().getNodeData().getReference().setDoi( value );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    formatError( mtn, e );
                    break;
                }
                break;
            case EVENTS_DUPLICATIONS:
                if ( !getMyNode().getNodeData().isHasEvent() ) {
                    getMyNode().getNodeData().setEvent( new Event() );
                }
                getMyNode().getNodeData().getEvent().setDuplications( parsePositiveInt( mtn, value ) );
                break;
            case EVENTS_SPECIATIONS:
                if ( !getMyNode().getNodeData().isHasEvent() ) {
                    getMyNode().getNodeData().setEvent( new Event() );
                }
                getMyNode().getNodeData().getEvent().setSpeciations( parsePositiveInt( mtn, value ) );
                break;
            case EVENTS_GENE_LOSSES:
                if ( !getMyNode().getNodeData().isHasEvent() ) {
                    getMyNode().getNodeData().setEvent( new Event() );
                }
                getMyNode().getNodeData().getEvent().setGeneLosses( parsePositiveInt( mtn, value ) );
                break;
            case DATE_DESCRIPTION:
                ForesterUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setDesc( value );
                break;
            case DATE_MAX:
                ForesterUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setMax( parseBigDecimal( mtn, value ) );
                break;
            case DATE_MIN:
                ForesterUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setMin( parseBigDecimal( mtn, value ) );
                break;
            case DATE_UNIT:
                ForesterUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setUnit( value );
                break;
            case DATE_VALUE:
                ForesterUtil.ensurePresenceOfDate( getMyNode() );
                getMyNode().getNodeData().getDate().setValue( parseBigDecimal( mtn, value ) );
                break;
            case DIST_ALT:
                ForesterUtil.ensurePresenceOfDistribution( getMyNode() );
                getMyNode().getNodeData().getDistribution().setAltitude( parseBigDecimal( mtn, value ) );
                break;
            case DIST_DESC:
                ForesterUtil.ensurePresenceOfDistribution( getMyNode() );
                getMyNode().getNodeData().getDistribution().setDescription( value );
                break;
            case DIST_GEODETIC:
                ForesterUtil.ensurePresenceOfDistribution( getMyNode() );
                getMyNode().getNodeData().getDistribution().setGeodeticDatum( value );
                break;
            case DIST_LAT:
                ForesterUtil.ensurePresenceOfDistribution( getMyNode() );
                getMyNode().getNodeData().getDistribution().setLatitude( parseBigDecimal( mtn, value ) );
                break;
            case DIST_LONG:
                ForesterUtil.ensurePresenceOfDistribution( getMyNode() );
                getMyNode().getNodeData().getDistribution().setLongitude( parseBigDecimal( mtn, value ) );
                break;
            default:
                throw new IllegalArgumentException( "unknown: " + tag );
        }
        getJTree().repaint();
        getTreePanel().setEdited( true );
        getTreePanel().repaint();
    }

    private enum PHYLOXML_TAG {
        NODE_NAME,
        NODE_BRANCH_LENGTH,
        TAXONOMY_CODE,
        TAXONOMY_SCIENTIFIC_NAME,
        TAXONOMY_AUTHORITY,
        TAXONOMY_COMMON_NAME,
        TAXONOMY_SYNONYM,
        TAXONOMY_RANK,
        TAXONOMY_URI,
        SEQ_SYMBOL,
        SEQ_NAME,
        SEQ_LOCATION,
        SEQ_TYPE,
        SEQ_MOL_SEQ,
        SEQ_URI,
        DATE_DESCRIPTION,
        DATE_VALUE,
        DATE_MIN,
        DATE_MAX,
        DATE_UNIT,
        TAXONOMY_ID_VALUE,
        TAXONOMY_ID_PROVIDER,
        SEQ_ACC_VALUE,
        SEQ_ACC_SOURCE,
        CONFIDENCE_VALUE,
        CONFIDENCE_TYPE,
        LIT_REFERENCE_DESC,
        LIT_REFERENCE_DOI,
        EVENTS_DUPLICATIONS,
        EVENTS_SPECIATIONS,
        EVENTS_GENE_LOSSES,
        DIST_DESC,
        DIST_GEODETIC,
        DIST_LAT,
        DIST_LONG,
        DIST_ALT
    }

    private class TagNumber {

        final private PHYLOXML_TAG _tag;
        final private int          _number;

        TagNumber( final PHYLOXML_TAG tag, final int number ) {
            _tag = tag;
            _number = number;
        }

        int getNumber() {
            return _number;
        }

        PHYLOXML_TAG getTag() {
            return _tag;
        }

        @Override
        public String toString() {
            return getTag() + "_" + getNumber();
        }
    }
}
