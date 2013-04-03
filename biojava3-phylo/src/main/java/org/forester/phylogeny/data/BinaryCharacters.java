// $Id: BinaryCharacters.java,v 1.21 2009/10/26 23:29:39 cmzmasek Exp $
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public class BinaryCharacters implements PhylogenyData {

    public final static int         COUNT_DEFAULT = -1;
    private final SortedSet<String> _present;
    private final SortedSet<String> _gained;
    private final SortedSet<String> _lost;
    private final int               _present_count;
    private final int               _gained_count;
    private final int               _lost_count;
    private String                  _type;

    public BinaryCharacters() {
        _present = new TreeSet<String>();
        _gained = new TreeSet<String>();
        _lost = new TreeSet<String>();
        _present_count = COUNT_DEFAULT;
        _gained_count = COUNT_DEFAULT;
        _lost_count = COUNT_DEFAULT;
    }

    public BinaryCharacters( final SortedSet<String> present_characters,
                             final SortedSet<String> gained_characters,
                             final SortedSet<String> lost_characters,
                             final String type ) {
        _present = present_characters;
        _gained = gained_characters;
        _lost = lost_characters;
        _type = type;
        _present_count = COUNT_DEFAULT;
        _gained_count = COUNT_DEFAULT;
        _lost_count = COUNT_DEFAULT;
    }

    public BinaryCharacters( final SortedSet<String> present_characters,
                             final SortedSet<String> gained_characters,
                             final SortedSet<String> lost_characters,
                             final String type,
                             final int present_count,
                             final int gained_count,
                             final int lost_count ) {
        _present = present_characters;
        _gained = gained_characters;
        _lost = lost_characters;
        _type = type;
        _present_count = present_count;
        _gained_count = gained_count;
        _lost_count = lost_count;
        validate();
    }

    private void addCharacters( final String indentation, final Writer w, final String[] present ) throws IOException {
        for( final String string : present ) {
            PhylogenyDataUtil.appendElement( w, PhyloXmlMapping.BINARY_CHARACTER, string, indentation );
        }
    }

    public void addGainedCharacter( final String binary_character ) {
        if ( getLostCharacters().contains( binary_character ) ) {
            throw new IllegalArgumentException( "attempt to add binary character [" + binary_character
                    + "] to gained characters but is already listed as lost" );
        }
        getGainedCharacters().add( binary_character );
    }

    public void addLostCharacter( final String binary_character ) {
        if ( getPresentCharacters().contains( binary_character ) ) {
            throw new IllegalArgumentException( "attempt to add binary character [" + binary_character
                    + "] to lost characters but is already listed as present" );
        }
        if ( getGainedCharacters().contains( binary_character ) ) {
            throw new IllegalArgumentException( "attempt to add binary character [" + binary_character
                    + "] to lost characters but is already listed as gained" );
        }
        getLostCharacters().add( binary_character );
    }

    public void addPresentCharacter( final String binary_character ) {
        if ( getLostCharacters().contains( binary_character ) ) {
            throw new IllegalArgumentException( "attempt to add binary character [" + binary_character
                    + "] to present characters but is already listed as lost" );
        }
        getPresentCharacters().add( binary_character );
    }

    @Override
    public StringBuffer asSimpleText() {
        return asText();
    }

    @Override
    public StringBuffer asText() {
        validate();
        final StringBuffer sb = new StringBuffer();
        sb.append( "present [" );
        sb.append( getPresentCount() );
        sb.append( "]: " );
        sb.append( getPresentCharactersAsStringBuffer() );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        sb.append( "gained  [ " );
        sb.append( getGainedCount() );
        sb.append( "]: " );
        sb.append( getGainedCharactersAsStringBuffer() );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        sb.append( "lost    [" );
        sb.append( getLostCount() );
        sb.append( "]: " );
        sb.append( getLostCharactersAsStringBuffer() );
        return sb;
    }

    @Override
    /**
     * Not a deep copy.
     * 
     */
    public PhylogenyData copy() {
        validate();
        return new BinaryCharacters( getPresentCharacters(),
                                     getGainedCharacters(),
                                     getLostCharacters(),
                                     getType(),
                                     getPresentCount(),
                                     getGainedCount(),
                                     getLostCount() );
    }

    public SortedSet<String> getGainedCharacters() {
        return _gained;
    }

    public String[] getGainedCharactersAsStringArray() {
        return sortedSetToStringArray( getGainedCharacters() );
    }

    public StringBuffer getGainedCharactersAsStringBuffer() {
        return sortedSetToStringBuffer( getGainedCharacters(), " " );
    }

    public int getGainedCount() {
        return _gained_count;
    }

    public SortedSet<String> getLostCharacters() {
        return _lost;
    }

    public String[] getLostCharactersAsStringArray() {
        return sortedSetToStringArray( getLostCharacters() );
    }

    public StringBuffer getLostCharactersAsStringBuffer() {
        return sortedSetToStringBuffer( getLostCharacters(), " " );
    }

    public int getLostCount() {
        return _lost_count;
    }

    public SortedSet<String> getPresentCharacters() {
        return _present;
    }

    public String[] getPresentCharactersAsStringArray() {
        return sortedSetToStringArray( getPresentCharacters() );
    }

    public StringBuffer getPresentCharactersAsStringBuffer() {
        return sortedSetToStringBuffer( getPresentCharacters(), " " );
    }

    public int getPresentCount() {
        return _present_count;
    }

    public String getType() {
        return _type;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    public void setType( final String type ) {
        _type = type;
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        validate();
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer,
                                      PhyloXmlMapping.BINARY_CHARACTERS,
                                      PhyloXmlMapping.BINARY_CHARACTERS_TYPE_ATTR,
                                      getType(),
                                      PhyloXmlMapping.BINARY_CHARACTERS_GAINED_COUNT_ATTR,
                                      getGainedCount() != COUNT_DEFAULT ? String.valueOf( getGainedCount() ) : "",
                                      PhyloXmlMapping.BINARY_CHARACTERS_LOST_COUNT_ATTR,
                                      getLostCount() != COUNT_DEFAULT ? String.valueOf( getLostCount() ) : "",
                                      PhyloXmlMapping.BINARY_CHARACTERS_PRESENT_COUNT_ATTR,
                                      getPresentCount() != COUNT_DEFAULT ? String.valueOf( getPresentCount() ) : "" );
        final String my_ind = indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE;
        if ( getGainedCharacters().size() > 0 ) {
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( my_ind );
            PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.BINARY_CHARACTERS_GAINED );
            addCharacters( my_ind, writer, getGainedCharactersAsStringArray() );
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( my_ind );
            PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.BINARY_CHARACTERS_GAINED );
        }
        if ( getLostCharacters().size() > 0 ) {
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( my_ind );
            PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.BINARY_CHARACTERS_LOST );
            addCharacters( my_ind, writer, getLostCharactersAsStringArray() );
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( my_ind );
            PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.BINARY_CHARACTERS_LOST );
        }
        if ( getPresentCharacters().size() > 0 ) {
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( my_ind );
            PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.BINARY_CHARACTERS_PRESENT );
            addCharacters( my_ind, writer, getPresentCharactersAsStringArray() );
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( my_ind );
            PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.BINARY_CHARACTERS_PRESENT );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.BINARY_CHARACTERS );
    }

    @Override
    public String toString() {
        return asText().toString();
    }

    private void validate() {
        if ( ( getPresentCount() != COUNT_DEFAULT ) && ( getPresentCharacters().size() > 0 )
                && ( getPresentCount() != getPresentCharacters().size() ) ) {
            throw new IllegalStateException( "present characters size and count are unequal" );
        }
        if ( ( getGainedCount() != COUNT_DEFAULT ) && ( getGainedCharacters().size() > 0 )
                && ( getGainedCount() != getGainedCharacters().size() ) ) {
            throw new IllegalStateException( "gained characters size and count are unequal" );
        }
        if ( ( getLostCount() != COUNT_DEFAULT ) && ( getLostCharacters().size() > 0 )
                && ( getLostCount() != getLostCharacters().size() ) ) {
            throw new IllegalStateException( "lost characters size and count are unequal" );
        }
    }

    private static String[] sortedSetToStringArray( final SortedSet<String> set ) {
        final String[] chars = new String[ set.size() ];
        final Iterator<String> it = set.iterator();
        int i = 0;
        while ( it.hasNext() ) {
            chars[ i++ ] = it.next();
        }
        return chars;
    }

    private static StringBuffer sortedSetToStringBuffer( final SortedSet<String> set, final String separator ) {
        final StringBuffer sb = new StringBuffer();
        final Iterator<String> it = set.iterator();
        while ( it.hasNext() ) {
            sb.append( it.next() );
            if ( it.hasNext() ) {
                sb.append( separator );
            }
        }
        return sb;
    }
}
