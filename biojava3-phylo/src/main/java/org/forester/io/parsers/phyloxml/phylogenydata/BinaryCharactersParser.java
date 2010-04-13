// $Id: BinaryCharactersParser.java,v 1.2 2009/11/03 19:16:34 cmzmasek Exp $
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

package org.forester.io.parsers.phyloxml.phylogenydata;

import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.data.PhylogenyData;

public class BinaryCharactersParser implements PhylogenyDataPhyloXmlParser {

    private static final BinaryCharactersParser _instance;
    static {
        try {
            _instance = new BinaryCharactersParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private BinaryCharactersParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhylogenyParserException {
        final SortedSet<String> present = new TreeSet<String>();
        final SortedSet<String> gained = new TreeSet<String>();
        final SortedSet<String> lost = new TreeSet<String>();
        String type = "";
        int present_count = BinaryCharacters.COUNT_DEFAULT;
        int gained_count = BinaryCharacters.COUNT_DEFAULT;
        int lost_count = BinaryCharacters.COUNT_DEFAULT;
        if ( element.isHasAttribute( PhyloXmlMapping.BINARY_CHARACTERS_TYPE_ATTR ) ) {
            type = element.getAttribute( PhyloXmlMapping.BINARY_CHARACTERS_TYPE_ATTR );
        }
        try {
            if ( element.isHasAttribute( PhyloXmlMapping.BINARY_CHARACTERS_PRESENT_COUNT_ATTR ) ) {
                present_count = Integer.parseInt( element
                        .getAttribute( PhyloXmlMapping.BINARY_CHARACTERS_PRESENT_COUNT_ATTR ) );
            }
            if ( element.isHasAttribute( PhyloXmlMapping.BINARY_CHARACTERS_GAINED_COUNT_ATTR ) ) {
                gained_count = Integer.parseInt( element
                        .getAttribute( PhyloXmlMapping.BINARY_CHARACTERS_GAINED_COUNT_ATTR ) );
            }
            if ( element.isHasAttribute( PhyloXmlMapping.BINARY_CHARACTERS_LOST_COUNT_ATTR ) ) {
                lost_count = Integer
                        .parseInt( element.getAttribute( PhyloXmlMapping.BINARY_CHARACTERS_LOST_COUNT_ATTR ) );
            }
        }
        catch ( final NumberFormatException e ) {
            throw new PhylogenyParserException( "failed to parse integer from element " + element.getQualifiedName() );
        }
        for( int i = 0; i < element.getNumberOfChildElements(); ++i ) {
            final XmlElement child_element = element.getChildElement( i );
            if ( child_element.getQualifiedName().equals( PhyloXmlMapping.BINARY_CHARACTERS_PRESENT ) ) {
                parseCharacters( present, child_element );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.BINARY_CHARACTERS_GAINED ) ) {
                parseCharacters( gained, child_element );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.BINARY_CHARACTERS_LOST ) ) {
                parseCharacters( lost, child_element );
            }
        }
        BinaryCharacters bc = null;
        if ( present_count != BinaryCharacters.COUNT_DEFAULT ) {
            bc = new BinaryCharacters( present, gained, lost, type, present_count, gained_count, lost_count );
        }
        else {
            bc = new BinaryCharacters( present, gained, lost, type );
        }
        return bc;
    }

    private void parseCharacters( final SortedSet<String> present, final XmlElement child_element ) {
        for( int j = 0; j < child_element.getNumberOfChildElements(); ++j ) {
            final XmlElement child_child_element = child_element.getChildElement( j );
            if ( child_child_element.getQualifiedName().equals( PhyloXmlMapping.BINARY_CHARACTER ) ) {
                present.add( child_child_element.getValueAsString() );
            }
        }
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
