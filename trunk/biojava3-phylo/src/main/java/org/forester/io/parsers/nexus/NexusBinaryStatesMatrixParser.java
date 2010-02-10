// $Id: NexusBinaryStatesMatrixParser.java,v 1.1 2009/02/04 01:10:14 cmzmasek
// Exp $
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2009-2010 Christian M. Zmasek
// Copyright (C) 2009-2010 Burnham Institute for Medical Research
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
// WWW: www.phylosoft.org/forester/

package org.forester.io.parsers.nexus;

import java.io.BufferedReader;
import java.io.IOException;

import org.forester.io.parsers.ParserUtils;
import org.forester.io.parsers.PhylogenyParserException;
import org.forester.phylogenyinference.BasicCharacterStateMatrix;
import org.forester.phylogenyinference.CharacterStateMatrix;
import org.forester.phylogenyinference.CharacterStateMatrix.BinaryStates;

public class NexusBinaryStatesMatrixParser {

    private Object                             _nexus_source;
    private CharacterStateMatrix<BinaryStates> _matrix;
    private int                                _nchar;
    private int                                _ntax;

    public CharacterStateMatrix<BinaryStates> getMatrix() {
        return _matrix;
    }

    public int getNChar() {
        return _nchar;
    }

    private Object getNexusSource() {
        return _nexus_source;
    }

    public int getNTax() {
        return _ntax;
    }

    public void parse() throws IOException {
        reset();
        final BufferedReader reader = ParserUtils.createReader( getNexusSource() );
        String line;
        boolean in_matrix = false;
        int identifier_index = 0;
        int max_character_index = -1;
        while ( ( line = reader.readLine() ) != null ) {
            line = line.trim();
            if ( ( line.length() > 0 ) && !line.startsWith( "#" ) && !line.startsWith( ">" ) ) {
                if ( line.toLowerCase().indexOf( NexusConstants.NCHAR.toLowerCase() ) >= 0 ) {
                    final int i = line.toLowerCase().indexOf( NexusConstants.NCHAR.toLowerCase() );
                    String s = line.toLowerCase().substring( i + 6 );
                    s = s.replace( ';', ' ' ).trim();
                    setNChar( Integer.parseInt( s ) );
                }
                else if ( line.toLowerCase().indexOf( NexusConstants.NTAX.toLowerCase() ) >= 0 ) {
                    final int i = line.toLowerCase().indexOf( NexusConstants.NTAX.toLowerCase() );
                    String s = line.toLowerCase().substring( i + 5 );
                    s = s.replace( ';', ' ' ).trim();
                    setNTax( Integer.parseInt( s ) );
                }
                else if ( line.toLowerCase().startsWith( NexusConstants.MATRIX.toLowerCase() ) ) {
                    in_matrix = true;
                    if ( getNTax() < 1 ) {
                        throw new NexusFormatException( "did not encounter " + NexusConstants.NTAX );
                    }
                    if ( getNChar() < 1 ) {
                        throw new NexusFormatException( "did not encounter " + NexusConstants.NCHAR );
                    }
                    if ( getMatrix() != null ) {
                        throw new NexusFormatException( "more than one matrix present" );
                    }
                    setMatrix( new BasicCharacterStateMatrix<BinaryStates>( getNTax(), getNChar() ) );
                }
                else if ( line.toLowerCase().startsWith( NexusConstants.END.toLowerCase() ) ) {
                    in_matrix = false;
                }
                else if ( in_matrix ) {
                    final String[] line_ary = line.split( "\\s+" );
                    final String label = line_ary[ 0 ].trim();
                    String states_str = line_ary[ 1 ].trim();
                    if ( states_str.endsWith( ";" ) ) {
                        in_matrix = false;
                        states_str = states_str.substring( 0, states_str.length() - 1 );
                    }
                    final char[] states = states_str.toCharArray();
                    getMatrix().setIdentifier( identifier_index, label );
                    int character_index = 0;
                    for( final char state : states ) {
                        if ( state == BinaryStates.PRESENT.toChar() ) {
                            try {
                                getMatrix().setState( identifier_index, character_index, BinaryStates.PRESENT );
                            }
                            catch ( final ArrayIndexOutOfBoundsException ex ) {
                                throw new NexusFormatException( "problem at line " + line + " [" + ex + "]" );
                            }
                        }
                        else if ( state == BinaryStates.ABSENT.toChar() ) {
                            try {
                                getMatrix().setState( identifier_index, character_index, BinaryStates.ABSENT );
                            }
                            catch ( final ArrayIndexOutOfBoundsException ex ) {
                                throw new NexusFormatException( "problem at line " + line + " [" + ex + "]" );
                            }
                        }
                        else {
                            throw new NexusFormatException( "illegal state " + state );
                        }
                        ++character_index;
                    }
                    if ( ( max_character_index > 0 ) && ( max_character_index != character_index ) ) {
                        throw new NexusFormatException( "unequal number of characters at line " + line );
                    }
                    max_character_index = character_index;
                    ++identifier_index;
                }
            }
        }
    }

    private void reset() {
        setMatrix( null );
        setNChar( -1 );
        setNTax( -1 );
    }

    private void setMatrix( final CharacterStateMatrix<BinaryStates> matrix ) {
        _matrix = matrix;
    }

    private void setNChar( final int nchar ) {
        _nchar = nchar;
    }

    private void setNTax( final int ntax ) {
        _ntax = ntax;
    }

    public void setSource( final Object nexus_source ) throws PhylogenyParserException, IOException {
        if ( nexus_source == null ) {
            throw new PhylogenyParserException( getClass() + ": attempt to parse null object." );
        }
        _nexus_source = nexus_source;
    }
}
