// $Id: TolXmlHandler.java,v 1.9 2009/11/17 19:53:23 cmzmasek Exp $
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

package org.forester.io.parsers.tol;

import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

public final class TolXmlHandler extends DefaultHandler {

    private String                    _current_element_name;
    private Phylogeny                 _current_phylogeny;
    private List<Phylogeny>           _phylogenies;
    private XmlElement                _current_xml_element;
    private PhylogenyNode             _current_node;
    private final static StringBuffer _buffer = new StringBuffer();

    TolXmlHandler() {
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
                && !getCurrentElementName().equals( TolXmlMapping.CLADE )
                && !getCurrentElementName().equals( TolXmlMapping.PHYLOGENY ) ) {
            getCurrentXmlElement().setValue( new String( chars, start_index, end_index ).trim() );
        }
    }

    @Override
    public void endElement( final String namespace_uri, final String local_name, final String qualified_name )
            throws SAXException {
        if ( ForesterUtil.isEmpty( namespace_uri ) || namespace_uri.startsWith( ForesterConstants.PHYLO_XML_LOCATION ) ) {
            if ( local_name.equals( TolXmlMapping.CLADE ) ) {
                try {
                    TolXmlHandler.mapElementToPhylogenyNode( getCurrentXmlElement(), getCurrentNode() );
                    if ( !getCurrentNode().isRoot() ) {
                        setCurrentNode( getCurrentNode().getParent() );
                    }
                    setCurrentXmlElement( getCurrentXmlElement().getParent() );
                }
                catch ( final PhylogenyParserException ex ) {
                    throw new SAXException( ex.getMessage() );
                }
            }
            else if ( local_name.equals( TolXmlMapping.PHYLOGENY ) ) {
                try {
                    TolXmlHandler.mapElementToPhylogeny( getCurrentXmlElement(), getCurrentPhylogeny() );
                }
                catch ( final PhylogenyParserException ex ) {
                    throw new SAXException( ex.getMessage() );
                }
                finishPhylogeny();
                reset();
            }
            else if ( ( getCurrentPhylogeny() != null ) && ( getCurrentXmlElement().getParent() != null ) ) {
                setCurrentXmlElement( getCurrentXmlElement().getParent() );
            }
            setCurrentElementName( null );
        }
    }

    private void finishPhylogeny() throws SAXException {
        getCurrentPhylogeny().setRooted( true );
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
        setCurrentElementName( local_name );
        if ( local_name.equals( TolXmlMapping.CLADE ) ) {
            final XmlElement element = new XmlElement( namespace_uri, local_name, local_name, attributes );
            getCurrentXmlElement().addChildElement( element );
            setCurrentXmlElement( element );
            newClade();
        }
        else if ( local_name.equals( TolXmlMapping.PHYLOGENY ) ) {
            setCurrentXmlElement( new XmlElement( "", "", "", null ) );
            newPhylogeny();
        }
        else if ( getCurrentPhylogeny() != null ) {
            final XmlElement element = new XmlElement( namespace_uri, local_name, local_name, attributes );
            getCurrentXmlElement().addChildElement( element );
            setCurrentXmlElement( element );
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
        // Not needed for now.
    }

    private static void mapElementToPhylogenyNode( final XmlElement xml_element, final PhylogenyNode node )
            throws PhylogenyParserException {
        if ( xml_element.isHasAttribute( TolXmlMapping.NODE_ID_ATTR ) ) {
            final String id = xml_element.getAttribute( TolXmlMapping.NODE_ID_ATTR );
            if ( !ForesterUtil.isEmpty( id ) ) {
                if ( !node.getNodeData().isHasTaxonomy() ) {
                    node.getNodeData().setTaxonomy( new Taxonomy() );
                }
                node.getNodeData().getTaxonomy()
                        .setIdentifier( new Identifier( id, TolXmlMapping.TOL_TAXONOMY_ID_TYPE ) );
            }
        }
        boolean put_into_scientific_name = false;
        if ( xml_element.isHasAttribute( TolXmlMapping.NODE_ITALICIZENAME_ATTR ) ) {
            final String ital = xml_element.getAttribute( TolXmlMapping.NODE_ITALICIZENAME_ATTR );
            if ( !ForesterUtil.isEmpty( ital ) && ital.equals( "1" ) ) {
                put_into_scientific_name = true;
            }
        }
        for( int i = 0; i < xml_element.getNumberOfChildElements(); ++i ) {
            final XmlElement element = xml_element.getChildElement( i );
            final String qualified_name = element.getQualifiedName();
            if ( qualified_name.equals( TolXmlMapping.TAXONOMY_NAME ) ) {
                final String name = element.getValueAsString();
                if ( !ForesterUtil.isEmpty( name ) ) {
                    if ( !node.getNodeData().isHasTaxonomy() ) {
                        node.getNodeData().setTaxonomy( new Taxonomy() );
                    }
                    if ( put_into_scientific_name ) {
                        node.getNodeData().getTaxonomy().setScientificName( name );
                    }
                    else {
                        node.getNodeData().getTaxonomy().setCommonName( name );
                    }
                }
            }
            else if ( qualified_name.equals( TolXmlMapping.AUTHORITY ) ) {
                String auth = element.getValueAsString();
                if ( !ForesterUtil.isEmpty( auth ) && !auth.equalsIgnoreCase( "null" ) ) {
                    if ( !node.getNodeData().isHasTaxonomy() ) {
                        node.getNodeData().setTaxonomy( new Taxonomy() );
                    }
                    auth = auth.replaceAll( "&amp;", "&" );
                    node.getNodeData().getTaxonomy().setAuthority( auth );
                }
            }
            else if ( qualified_name.equals( TolXmlMapping.AUTHDATE ) ) {
                final String authdate = element.getValueAsString();
                if ( !ForesterUtil.isEmpty( authdate ) && !authdate.equalsIgnoreCase( "null" ) ) {
                    if ( node.getNodeData().isHasTaxonomy()
                            && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getAuthority() ) ) {
                        _buffer.setLength( 0 );
                        _buffer.append( node.getNodeData().getTaxonomy().getAuthority() );
                        _buffer.append( " " );
                        _buffer.append( authdate );
                        node.getNodeData().getTaxonomy().setAuthority( _buffer.toString() );
                    }
                }
            }
            else if ( qualified_name.equals( TolXmlMapping.OTHERNAMES ) ) {
                for( int j = 0; j < element.getNumberOfChildElements(); ++j ) {
                    final XmlElement element_j = element.getChildElement( j );
                    if ( element_j.getQualifiedName().equals( TolXmlMapping.OTHERNAME ) ) {
                        for( int z = 0; z < element_j.getNumberOfChildElements(); ++z ) {
                            final XmlElement element_z = element_j.getChildElement( z );
                            if ( element_z.getQualifiedName().equals( TolXmlMapping.OTHERNAME_NAME ) ) {
                                final String syn = element_z.getValueAsString();
                                if ( !ForesterUtil.isEmpty( syn ) && !syn.equalsIgnoreCase( "null" ) ) {
                                    if ( !node.getNodeData().isHasTaxonomy() ) {
                                        node.getNodeData().setTaxonomy( new Taxonomy() );
                                    }
                                    node.getNodeData().getTaxonomy().getSynonyms().add( syn );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}