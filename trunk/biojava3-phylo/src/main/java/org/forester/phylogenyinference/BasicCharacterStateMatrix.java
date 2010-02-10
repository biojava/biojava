// $Id: BasicCharacterStateMatrix.java,v 1.33 2009/10/26 23:29:39 cmzmasek Exp $
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

package org.forester.phylogenyinference;

import java.io.IOException;
import java.io.Writer;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.forester.io.parsers.nexus.NexusConstants;
import org.forester.util.ForesterUtil;

public class BasicCharacterStateMatrix<S> implements CharacterStateMatrix<S> {

    final Object[][]           _states;
    final String[]             _identifiers;
    final String[]             _characters;
    final Map<String, Integer> _identifier_index_map;
    final Map<String, Integer> _character_index_map;

    public BasicCharacterStateMatrix( final int number_of_identifiers, final int number_of_characters ) {
        _states = new Object[ number_of_identifiers ][ number_of_characters ];
        _identifiers = new String[ number_of_identifiers ];
        _characters = new String[ number_of_characters ];
        _identifier_index_map = new HashMap<String, Integer>( number_of_identifiers );
        _character_index_map = new HashMap<String, Integer>( number_of_characters );
    }

    public BasicCharacterStateMatrix( final int number_of_identifiers,
                                      final int number_of_characters,
                                      final S default_state ) {
        this( number_of_identifiers, number_of_identifiers );
        for( int identifier = 0; identifier < number_of_identifiers; ++identifier ) {
            for( int character = 0; character < number_of_characters; ++character ) {
                setState( identifier, character, default_state );
            }
        }
    }

    public BasicCharacterStateMatrix( final List<List<S>> states ) {
        if ( ( states == null ) || ( states.size() < 1 ) || ( states.get( 0 ) == null ) ) {
            throw new IllegalArgumentException( "attempt to create character state matrix from empty list" );
        }
        final int number_of_characters = states.get( 0 ).size();
        final int number_of_identifiers = states.size();
        _states = new Object[ number_of_identifiers ][ number_of_characters ];
        _identifiers = new String[ number_of_identifiers ];
        _characters = new String[ number_of_characters ];
        _identifier_index_map = new HashMap<String, Integer>( number_of_identifiers );
        _character_index_map = new HashMap<String, Integer>( number_of_characters );
        for( int identifier = 0; identifier < number_of_identifiers; ++identifier ) {
            for( int character = 0; character < number_of_characters; ++character ) {
                setState( identifier, character, states.get( identifier ).get( character ) );
            }
        }
    }

    public BasicCharacterStateMatrix( final S[][] states ) {
        this( states.length, states[ 0 ].length );
        for( int identifier = 0; identifier < states.length; ++identifier ) {
            for( int character = 0; character < states[ 0 ].length; ++character ) {
                setState( identifier, character, states[ identifier ][ character ] );
            }
        }
    }

    @Override
    public boolean containsCharacter( final String character ) {
        return _character_index_map.containsKey( character );
    }

    @Override
    public boolean containsIdentifier( final String identifier ) {
        return _identifier_index_map.containsKey( identifier );
    }

    public CharacterStateMatrix<S> copy() {
        final CharacterStateMatrix<S> new_matrix = new BasicCharacterStateMatrix<S>( getNumberOfIdentifiers(),
                                                                                     getNumberOfCharacters() );
        for( int character = 0; character < getNumberOfCharacters(); ++character ) {
            if ( getCharacter( character ) != null ) {
                new_matrix.setCharacter( character, getCharacter( character ) );
            }
        }
        for( int identifier = 0; identifier < getNumberOfIdentifiers(); ++identifier ) {
            if ( getIdentifier( identifier ) != null ) {
                new_matrix.setIdentifier( identifier, getIdentifier( identifier ) );
            }
            for( int character = 0; character < getNumberOfCharacters(); ++character ) {
                new_matrix.setState( identifier, character, getState( identifier, character ) );
            }
        }
        return new_matrix;
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check character state matrix equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check character state matrix to " + o + " [" + o.getClass()
                    + "]" );
        }
        else {
            final CharacterStateMatrix<S> other = ( CharacterStateMatrix<S> ) o;
            if ( ( getNumberOfIdentifiers() != other.getNumberOfIdentifiers() )
                    || ( getNumberOfCharacters() != other.getNumberOfCharacters() ) ) {
            }
            for( int identifier = 0; identifier < getNumberOfIdentifiers(); ++identifier ) {
                for( int character = 0; character < getNumberOfCharacters(); ++character ) {
                    final S s = getState( identifier, character );
                    final S os = other.getState( identifier, character );
                    if ( s == os ) {
                        continue;
                    }
                    else if ( ( s == null ) && ( os != null ) ) {
                        return false;
                    }
                    else if ( ( s != null ) && ( os == null ) ) {
                        return false;
                    }
                    else if ( !s.equals( other.getState( identifier, character ) ) ) {
                        return false;
                    }
                }
            }
            return true;
        }
    }

    @Override
    public String getCharacter( final int character_index ) {
        return _characters[ character_index ];
    }

    @Override
    public int getCharacterIndex( final String character ) {
        if ( !_character_index_map.containsKey( character ) ) {
            throw new IllegalArgumentException( "character [" + character + "] not found" );
        }
        return _character_index_map.get( character );
    }

    @Override
    public String getIdentifier( final int identifier_index ) {
        return _identifiers[ identifier_index ];
    }

    @Override
    public int getIdentifierIndex( final String identifier ) {
        if ( !_identifier_index_map.containsKey( identifier ) ) {
            throw new IllegalArgumentException( "indentifier [" + identifier + "] not found" );
        }
        return _identifier_index_map.get( identifier );
    }

    private int getLengthOfLongestState() {
        int longest = 0;
        for( int identifier = 0; identifier < getNumberOfIdentifiers(); ++identifier ) {
            for( int character = 0; character < getNumberOfCharacters(); ++character ) {
                final S s = getState( identifier, character );
                if ( s != null ) {
                    final int l = getState( identifier, character ).toString().length();
                    if ( l > longest ) {
                        longest = l;
                    }
                }
            }
        }
        return longest;
    }

    @Override
    public int getNumberOfCharacters() {
        if ( !isEmpty() ) {
            return _states[ 0 ].length;
        }
        else {
            return 0;
        }
    }

    @Override
    public int getNumberOfIdentifiers() {
        return _states.length;
    }

    @Override
    public S getState( final int identifier_index, final int character_index ) {
        return ( S ) _states[ identifier_index ][ character_index ];
    }

    @Override
    public S getState( final String identifier, final int character_index ) {
        if ( !containsIdentifier( identifier ) ) {
            throw new IllegalArgumentException( "identifier [" + identifier + "] not found" );
        }
        return getState( _identifier_index_map.get( identifier ), character_index );
    }

    @Override
    public S getState( final String identifier, final String character ) {
        if ( !containsIdentifier( identifier ) ) {
            throw new IllegalArgumentException( "identifier [" + identifier + "] not found" );
        }
        if ( !containsCharacter( character ) ) {
            throw new IllegalArgumentException( "character [" + character + "] not found" );
        }
        return getState( _identifier_index_map.get( identifier ), _character_index_map.get( character ) );
    }

    @Override
    public boolean isEmpty() {
        return getNumberOfIdentifiers() <= 0;
    }

    public CharacterStateMatrix<S> pivot() {
        final CharacterStateMatrix<S> new_matrix = new BasicCharacterStateMatrix<S>( getNumberOfCharacters(),
                                                                                     getNumberOfIdentifiers() );
        for( int character = 0; character < getNumberOfCharacters(); ++character ) {
            if ( getCharacter( character ) != null ) {
                new_matrix.setIdentifier( character, getCharacter( character ) );
            }
        }
        for( int identifier = 0; identifier < getNumberOfIdentifiers(); ++identifier ) {
            if ( getIdentifier( identifier ) != null ) {
                new_matrix.setCharacter( identifier, getIdentifier( identifier ) );
            }
            for( int character = 0; character < getNumberOfCharacters(); ++character ) {
                new_matrix.setState( character, identifier, getState( identifier, character ) );
            }
        }
        return new_matrix;
    }

    @Override
    public void setCharacter( final int character_index, final String character ) {
        if ( character == null ) {
            throw new IllegalArgumentException( "attempt to use null character" );
        }
        _characters[ character_index ] = character;
        if ( _character_index_map.containsKey( character ) ) {
            throw new IllegalArgumentException( "character [" + character + "] is not unique" );
        }
        _character_index_map.put( character, character_index );
    }

    @Override
    public void setIdentifier( final int identifier_index, final String identifier ) {
        if ( identifier == null ) {
            throw new IllegalArgumentException( "attempt to use null identifier" );
        }
        _identifiers[ identifier_index ] = identifier;
        if ( _identifier_index_map.containsKey( identifier ) ) {
            throw new IllegalArgumentException( "identifier [" + identifier + "] is not unique" );
        }
        _identifier_index_map.put( identifier, identifier_index );
    }

    @Override
    public void setState( final int identifier_index, final int character_index, final S state ) {
        _states[ identifier_index ][ character_index ] = state;
    }

    @Override
    public void setState( final String identifier, final int character_index, final S state ) {
        if ( !_identifier_index_map.containsKey( identifier ) ) {
            throw new IllegalArgumentException( "identifier [" + identifier + "] not found" );
        }
        setState( _identifier_index_map.get( identifier ), character_index, state );
    }

    @Override
    public void setState( final String identifier, final String character, final S state ) {
        if ( !containsIdentifier( identifier ) ) {
            throw new IllegalArgumentException( "identifier [" + identifier + "] not found" );
        }
        if ( !containsCharacter( character ) ) {
            throw new IllegalArgumentException( "character [" + character + "] not found" );
        }
        setState( _identifier_index_map.get( identifier ), _character_index_map.get( character ), state );
    }

    private void toForester( final Writer writer ) throws IOException {
        final int longest = getLengthOfLongestState() + 5;
        writer.write( "Identifiers: " );
        writer.write( String.valueOf( getNumberOfIdentifiers() ) );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( "Characters : " );
        writer.write( String.valueOf( getNumberOfCharacters() ) );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( ForesterUtil.pad( "", 20, ' ', false ).toString() );
        writer.write( ' ' );
        for( int character = 0; character < getNumberOfCharacters(); ++character ) {
            final String c = getCharacter( character );
            writer.write( c != null ? ForesterUtil.pad( c, longest, ' ', false ).toString() : ForesterUtil
                    .pad( "", longest, ' ', false ).toString() );
            if ( character < getNumberOfCharacters() - 1 ) {
                writer.write( ' ' );
            }
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        for( int identifier = 0; identifier < getNumberOfIdentifiers(); ++identifier ) {
            if ( getIdentifier( identifier ) != null ) {
                writer.write( ForesterUtil.pad( getIdentifier( identifier ), 20, ' ', false ).toString() );
                writer.write( ' ' );
            }
            for( int character = 0; character < getNumberOfCharacters(); ++character ) {
                final S state = getState( identifier, character );
                writer.write( state != null ? ForesterUtil.pad( state.toString(), longest, ' ', false ).toString()
                        : ForesterUtil.pad( "", longest, ' ', false ).toString() );
                if ( character < getNumberOfCharacters() - 1 ) {
                    writer.write( ' ' );
                }
            }
            if ( identifier < getNumberOfIdentifiers() - 1 ) {
                writer.write( ForesterUtil.LINE_SEPARATOR );
            }
        }
    }

    private void toNexus( final Writer writer ) throws IOException {
        if ( isEmpty() ) {
            return;
        }
        writer.write( NexusConstants.NEXUS );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writeNexusTaxaBlock( writer );
        writeNexusBinaryChractersBlock( writer );
    }

    private void toPhylip( final Writer writer ) throws IOException {
        final int pad = 6;
        writer.write( ' ' );
        writer.write( ' ' );
        writer.write( ' ' );
        writer.write( ' ' );
        writer.write( getNumberOfIdentifiers() );
        writer.write( ' ' );
        writer.write( getNumberOfCharacters() );
        writer.write( ForesterUtil.LINE_SEPARATOR );
        for( int identifier = 0; identifier < getNumberOfIdentifiers(); ++identifier ) {
            if ( !ForesterUtil.isEmpty( getIdentifier( identifier ) ) ) {
                writer.write( ForesterUtil.pad( getIdentifier( identifier ), pad, ' ', false ).toString() );
                writer.write( ' ' );
                writer.write( ' ' );
            }
            else {
                throw new IllegalStateException( "PHYLIP format does not allow empty identifiers" );
            }
            writer.write( "" );
            for( int character = 0; character < getNumberOfCharacters(); ++character ) {
                final String state = getState( identifier, character ).toString();
                writer.write( state != null ? ForesterUtil.pad( state, pad, ' ', false ).toString() : ForesterUtil
                        .pad( "", pad, ' ', false ).toString() );
                if ( character < getNumberOfCharacters() - 1 ) {
                    writer.write( ' ' );
                    writer.write( ' ' );
                }
            }
            if ( identifier < getNumberOfIdentifiers() - 1 ) {
                writer.write( ForesterUtil.LINE_SEPARATOR );
            }
        }
    }

    //TODO
    //to format for microarray-style clustering
    // states are ints in this case
    //TODO
    public void toWriter( final Writer writer ) throws IOException {
        toForester( writer );
    }

    @Override
    public void toWriter( final Writer writer, final Format format ) throws IOException {
        switch ( format ) {
            case PHYLIP:
                toPhylip( writer );
                break;
            case FORESTER:
                toForester( writer );
                break;
            case NEXUS_BINARY:
                toNexus( writer );
                break;
            default:
                throw new IllegalArgumentException( "Unknown format:" + format );
        }
    }

    public void writeNexusBinaryChractersBlock( final Writer w ) throws IOException {
        //BEGIN CHARACTERS;
        // DIMENSIONS NCHAR=x;
        //BEGIN CHARSTATELABELS 
        // 1 bcl,
        // 2 tir,
        //END;
        // FORMAT DATATYPE=STANDARD SYMBOLS=;
        // MATRIX
        //  fish d d f
        //  frog s d f f
        //  snake x x x x;
        // END;
        w.write( NexusConstants.BEGIN_CHARACTERS );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( " " );
        w.write( NexusConstants.DIMENSIONS );
        w.write( " " );
        w.write( NexusConstants.NCHAR );
        w.write( "=" );
        w.write( String.valueOf( getNumberOfCharacters() ) );
        w.write( ";" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        writeNexusCharstatelabels( w );
        w.write( " " );
        w.write( NexusConstants.FORMAT );
        w.write( " " );
        w.write( NexusConstants.DATATYPE );
        w.write( "=" );
        w.write( NexusConstants.STANDARD );
        w.write( " " );
        w.write( NexusConstants.SYMBOLS );
        w.write( "=\"" );
        w.write( String.valueOf( BinaryStates.ABSENT ) );
        w.write( String.valueOf( BinaryStates.PRESENT ) );
        w.write( "\";" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        writeNexusMatrix( w );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( NexusConstants.END );
        w.write( ForesterUtil.LINE_SEPARATOR );
    }

    public void writeNexusCharstatelabels( final Writer w ) throws IOException {
        w.write( " " );
        w.write( NexusConstants.CHARSTATELABELS );
        w.write( ForesterUtil.LINE_SEPARATOR );
        for( int i = 0; i < getNumberOfCharacters(); ++i ) {
            w.write( "  " + ( i + 1 ) + " '" );
            w.write( getCharacter( i ) );
            w.write( "'" );
            if ( i < getNumberOfCharacters() - 1 ) {
                w.write( "," );
                w.write( ForesterUtil.LINE_SEPARATOR );
            }
        }
        w.write( ";" );
        w.write( ForesterUtil.LINE_SEPARATOR );
    }

    public void writeNexusMatrix( final Writer w ) throws IOException {
        w.write( " " );
        w.write( NexusConstants.MATRIX );
        w.write( ForesterUtil.LINE_SEPARATOR );
        for( int identifier = 0; identifier < getNumberOfIdentifiers(); ++identifier ) {
            if ( getIdentifier( identifier ) != null ) {
                w.write( "  " );
                w.write( ForesterUtil.pad( getIdentifier( identifier ), 20, ' ', false ).toString() );
                w.write( ' ' );
            }
            for( int character = 0; character < getNumberOfCharacters(); ++character ) {
                final S state = getState( identifier, character );
                if ( state == null ) {
                    throw new IllegalStateException( "character state matrix cannot contain null if to be represented in nexus format" );
                }
                if ( !( state instanceof BinaryStates ) ) {
                    throw new IllegalStateException( "nexus format representation expects binary character data - got ["
                            + getState( 0, 0 ).getClass() + "] instead" );
                }
                if ( state == BinaryStates.UNKNOWN ) {
                    throw new IllegalStateException( "character state matrix cannot contain unknown states if to be represented in nexus format" );
                }
                w.write( state.toString() );
            }
            if ( identifier < getNumberOfIdentifiers() - 1 ) {
                w.write( ForesterUtil.LINE_SEPARATOR );
            }
        }
        w.write( ";" );
    }

    public void writeNexusTaxaBlock( final Writer w ) throws IOException {
        //BEGIN TAXA;
        // DIMENSIONS NTAX=n;
        // TAXLABELS fish frog snake;
        //END;
        w.write( NexusConstants.BEGIN_TAXA );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( " " );
        w.write( NexusConstants.DIMENSIONS );
        w.write( " " );
        w.write( NexusConstants.NTAX );
        w.write( "=" );
        w.write( String.valueOf( getNumberOfIdentifiers() ) );
        w.write( ";" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( " " );
        w.write( NexusConstants.TAXLABELS );
        for( int i = 0; i < getNumberOfIdentifiers(); ++i ) {
            w.write( " " );
            w.write( getIdentifier( i ) );
        }
        w.write( ";" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( NexusConstants.END );
        w.write( ForesterUtil.LINE_SEPARATOR );
    }
}
