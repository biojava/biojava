// $Id: PhyloXmlHandler.java,v 1.6 2009/11/09 22:42:48 cmzmasek Exp $
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

package org.forester.io.parsers.phyloxml;

import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.phyloxml.phylogenydata.BinaryCharactersParser;
import org.forester.io.parsers.phyloxml.phylogenydata.BranchWidthParser;
import org.forester.io.parsers.phyloxml.phylogenydata.ColorParser;
import org.forester.io.parsers.phyloxml.phylogenydata.ConfidenceParser;
import org.forester.io.parsers.phyloxml.phylogenydata.DateParser;
import org.forester.io.parsers.phyloxml.phylogenydata.EventParser;
import org.forester.io.parsers.phyloxml.phylogenydata.IdentifierParser;
import org.forester.io.parsers.phyloxml.phylogenydata.PropertyParser;
import org.forester.io.parsers.phyloxml.phylogenydata.ReferenceParser;
import org.forester.io.parsers.phyloxml.phylogenydata.SequenceParser;
import org.forester.io.parsers.phyloxml.phylogenydata.TaxonomyParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PropertiesMap;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Reference;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

public final class PhyloXmlHandler extends DefaultHandler {

    private static final String PHYLOXML = "phyloxml";
    private String              _current_element_name;
    private Phylogeny           _current_phylogeny;
    private List<Phylogeny>     _phylogenies;
    private XmlElement          _current_xml_element;
    private PhylogenyNode       _current_node;

    PhyloXmlHandler() {
        // Constructor.
    }

    private void addNode() {
        final PhylogenyNode new_node = new PhylogenyNode();
        getCurrentNode().addAsChild( new_node );
        setCurrentNode( new_node );
    }

    @Override
    public void characters( final char[] chars, final int start_index, final int end_index ) {
        if ( ( ( getCurrentXmlElement() != null ) && ( getCurrentElementName() != null ) )
                && !getCurrentElementName().equals( PhyloXmlMapping.CLADE )
                && !getCurrentElementName().equals( PhyloXmlMapping.PHYLOGENY ) ) {
            if ( !ForesterUtil.isEmpty( getCurrentXmlElement().getValueAsString() ) ) {
                getCurrentXmlElement().appendValue( new String( chars, start_index, end_index ) );
            }
            else {
                getCurrentXmlElement().setValue( new String( chars, start_index, end_index ) );
            }
        }
    }

    @Override
    public void endElement( final String namespace_uri, final String local_name, final String qualified_name )
            throws SAXException {
        if ( ForesterUtil.isEmpty( namespace_uri ) || namespace_uri.startsWith( ForesterConstants.PHYLO_XML_LOCATION ) ) {
            if ( local_name.equals( PhyloXmlMapping.CLADE ) ) {
                try {
                    PhyloXmlHandler.mapElementToPhylogenyNode( getCurrentXmlElement(), getCurrentNode() );
                    if ( !getCurrentNode().isRoot() ) {
                        setCurrentNode( getCurrentNode().getParent() );
                    }
                    getCurrentXmlElement().setValue( null );
                    setCurrentXmlElement( getCurrentXmlElement().getParent() );
                }
                catch ( final PhylogenyParserException ex ) {
                    throw new SAXException( ex.getMessage() );
                }
            }
            else if ( local_name.equals( PhyloXmlMapping.PHYLOGENY ) ) {
                try {
                    PhyloXmlHandler.mapElementToPhylogeny( getCurrentXmlElement(), getCurrentPhylogeny() );
                }
                catch ( final PhylogenyParserException ex ) {
                    throw new SAXException( ex.getMessage() );
                }
                finishPhylogeny();
                reset();
            }
            else if ( local_name.equals( PHYLOXML ) ) {
                // Do nothing.
            }
            else if ( ( getCurrentPhylogeny() != null ) && ( getCurrentXmlElement().getParent() != null ) ) {
                setCurrentXmlElement( getCurrentXmlElement().getParent() );
            }
            setCurrentElementName( null );
        }
    }

    private void finishPhylogeny() throws SAXException {
        getCurrentPhylogeny().recalculateNumberOfExternalDescendants( false );
        getPhylogenies().add( getCurrentPhylogeny() );
    }

    private String getCurrentElementName() {
        return _current_element_name;
    }

    private PhylogenyNode getCurrentNode() {
        return _current_node;
    }

    private Phylogeny getCurrentPhylogeny() {
        return _current_phylogeny;
    }

    private XmlElement getCurrentXmlElement() {
        return _current_xml_element;
    }

    List<Phylogeny> getPhylogenies() {
        return _phylogenies;
    }

    private void init() {
        reset();
        setPhylogenies( new ArrayList<Phylogeny>() );
    }

    private void initCurrentNode() {
        if ( getCurrentNode() != null ) {
            throw new IllegalStateException( "attempt to create new current node when current node already exists" );
        }
        if ( getCurrentPhylogeny() == null ) {
            throw new IllegalStateException( "attempt to create new current node for non-existing phylogeny" );
        }
        final PhylogenyNode node = new PhylogenyNode();
        getCurrentPhylogeny().setRoot( node );
        setCurrentNode( getCurrentPhylogeny().getRoot() );
    }

    private void newClade() {
        if ( getCurrentNode() == null ) {
            initCurrentNode();
        }
        else {
            addNode();
        }
    }

    private void newPhylogeny() {
        setCurrentPhylogeny( new Phylogeny() );
    }

    private void reset() {
        setCurrentPhylogeny( null );
        setCurrentNode( null );
        setCurrentElementName( null );
        setCurrentXmlElement( null );
    }

    private void setCurrentElementName( final String element_name ) {
        _current_element_name = element_name;
    }

    private void setCurrentNode( final PhylogenyNode current_node ) {
        _current_node = current_node;
    }

    private void setCurrentPhylogeny( final Phylogeny phylogeny ) {
        _current_phylogeny = phylogeny;
    }

    private void setCurrentXmlElement( final XmlElement element ) {
        _current_xml_element = element;
    }

    private void setPhylogenies( final List<Phylogeny> phylogenies ) {
        _phylogenies = phylogenies;
    }

    @Override
    public void startDocument() throws SAXException {
        init();
    }

    @Override
    public void startElement( final String namespace_uri,
                              final String local_name,
                              final String qualified_name,
                              final Attributes attributes ) throws SAXException {
        if ( ForesterUtil.isEmpty( namespace_uri ) || namespace_uri.startsWith( ForesterConstants.PHYLO_XML_LOCATION ) ) {
            setCurrentElementName( local_name );
            if ( local_name.equals( PhyloXmlMapping.CLADE ) ) {
                final XmlElement element = new XmlElement( namespace_uri, local_name, local_name, attributes );
                getCurrentXmlElement().addChildElement( element );
                setCurrentXmlElement( element );
                newClade();
            }
            else if ( local_name.equals( PhyloXmlMapping.PHYLOGENY ) ) {
                setCurrentXmlElement( new XmlElement( "", "", "", null ) );
                newPhylogeny();
                final XmlElement element = new XmlElement( namespace_uri, local_name, local_name, attributes );
                if ( element.isHasAttribute( PhyloXmlMapping.PHYLOGENY_IS_REROOTABLE_ATTR ) ) {
                    getCurrentPhylogeny().setRerootable( Boolean.parseBoolean( element
                            .getAttribute( PhyloXmlMapping.PHYLOGENY_IS_REROOTABLE_ATTR ) ) );
                }
                if ( element.isHasAttribute( PhyloXmlMapping.PHYLOGENY_BRANCHLENGTH_UNIT_ATTR ) ) {
                    getCurrentPhylogeny().setDistanceUnit( element
                            .getAttribute( PhyloXmlMapping.PHYLOGENY_BRANCHLENGTH_UNIT_ATTR ) );
                }
                if ( element.isHasAttribute( PhyloXmlMapping.PHYLOGENY_IS_ROOTED_ATTR ) ) {
                    getCurrentPhylogeny().setRooted( Boolean.parseBoolean( element
                            .getAttribute( PhyloXmlMapping.PHYLOGENY_IS_ROOTED_ATTR ) ) );
                }
                if ( element.isHasAttribute( PhyloXmlMapping.PHYLOGENY_TYPE_ATTR ) ) {
                    getCurrentPhylogeny().setType( ( element.getAttribute( PhyloXmlMapping.PHYLOGENY_TYPE_ATTR ) ) );
                }
            }
            else if ( local_name.equals( PHYLOXML ) ) {
            }
            else if ( getCurrentPhylogeny() != null ) {
                final XmlElement element = new XmlElement( namespace_uri, local_name, local_name, attributes );
                getCurrentXmlElement().addChildElement( element );
                setCurrentXmlElement( element );
            }
        }
    }

    public static boolean attributeEqualsValue( final XmlElement element,
                                                final String attributeName,
                                                final String attributeValue ) {
        final String attr = element.getAttribute( attributeName );
        return ( ( attr != null ) && attr.equals( attributeValue ) );
    }

    public static String getAtttributeValue( final XmlElement element, final String attributeName ) {
        final String attr = element.getAttribute( attributeName );
        if ( attr != null ) {
            return attr;
        }
        else {
            return "";
        }
    }

    private static void mapElementToPhylogeny( final XmlElement xml_element, final Phylogeny phylogeny )
            throws PhylogenyParserException {
        for( int i = 0; i < xml_element.getNumberOfChildElements(); ++i ) {
            final XmlElement element = xml_element.getChildElement( i );
            final String qualified_name = element.getQualifiedName();
            if ( qualified_name.equals( PhyloXmlMapping.PHYLOGENY_NAME ) ) {
                phylogeny.setName( element.getValueAsString() );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.PHYLOGENY_DESCRIPTION ) ) {
                phylogeny.setDescription( element.getValueAsString() );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.IDENTIFIER ) ) {
                phylogeny.setIdentifier( ( Identifier ) IdentifierParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.CONFIDENCE ) ) {
                phylogeny.setConfidence( ( Confidence ) ConfidenceParser.getInstance().parse( element ) );
            }
        }
    }

    private static void mapElementToPhylogenyNode( final XmlElement xml_element, final PhylogenyNode node )
            throws PhylogenyParserException {
        if ( xml_element.isHasAttribute( PhyloXmlMapping.BRANCH_LENGTH ) ) {
            double d = 0;
            try {
                d = Double.parseDouble( xml_element.getAttribute( PhyloXmlMapping.BRANCH_LENGTH ) );
            }
            catch ( final NumberFormatException e ) {
                throw new PhylogenyParserException( "ill formatted distance in clade attribute ["
                        + xml_element.getAttribute( PhyloXmlMapping.BRANCH_LENGTH ) + "]: " + e.getMessage() );
            }
            node.setDistanceToParent( d );
        }
        for( int i = 0; i < xml_element.getNumberOfChildElements(); ++i ) {
            final XmlElement element = xml_element.getChildElement( i );
            final String qualified_name = element.getQualifiedName();
            if ( qualified_name.equals( PhyloXmlMapping.BRANCH_LENGTH ) ) {
                if ( node.getDistanceToParent() != PhylogenyNode.DISTANCE_DEFAULT ) {
                    throw new PhylogenyParserException( "ill advised attempt to set distance twice for the same clade (probably via element and via attribute)" );
                }
                node.setDistanceToParent( element.getValueAsDouble() );
            }
            if ( qualified_name.equals( PhyloXmlMapping.NODE_NAME ) ) {
                node.setName( element.getValueAsString() );
            }
            //  else if ( qualified_name.equals( PhyloXmlMapping.NODE_IDENTIFIER ) ) {
            //      node.getNodeData().setNodeIdentifier( ( Identifier ) IdentifierParser.getInstance().parse( element ) );
            //  }
            else if ( qualified_name.equals( PhyloXmlMapping.TAXONOMY ) ) {
                node.getNodeData().setTaxonomy( ( Taxonomy ) TaxonomyParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.SEQUENCE ) ) {
                node.getNodeData().setSequence( ( Sequence ) SequenceParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.DISTRIBUTION ) ) {
                //node.getNodeData().setDistribution( ( Distribution ) DistributionParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.CLADE_DATE ) ) {
                node.getNodeData().setDate( ( Date ) DateParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.REFERENCE ) ) {
                node.getNodeData().setReference( ( Reference ) ReferenceParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.BINARY_CHARACTERS ) ) {
                node.getNodeData().setBinaryCharacters( ( BinaryCharacters ) BinaryCharactersParser.getInstance()
                        .parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.COLOR ) ) {
                node.getBranchData().setBranchColor( ( BranchColor ) ColorParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.CONFIDENCE ) ) {
                node.getBranchData().addConfidence( ( Confidence ) ConfidenceParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.WIDTH ) ) {
                node.getBranchData().setBranchWidth( ( BranchWidth ) BranchWidthParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.EVENTS ) ) {
                node.getNodeData().setEvent( ( Event ) EventParser.getInstance().parse( element ) );
            }
            else if ( qualified_name.equals( PhyloXmlMapping.PROPERTY ) ) {
                if ( !node.getNodeData().isHasProperties() ) {
                    node.getNodeData().setProperties( new PropertiesMap() );
                }
                node.getNodeData().getProperties().addProperty( ( Property ) PropertyParser.getInstance()
                        .parse( element ) );
            }
        }
    }
}
