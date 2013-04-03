// $Id: PaupLogParser.java,v 1.3 2009/10/26 23:29:40 cmzmasek Exp $
//
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

package org.forester.io.parsers.nexus;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.ParserUtils;
import org.forester.io.parsers.PhylogenyParserException;
import org.forester.phylogenyinference.BasicCharacterStateMatrix;
import org.forester.phylogenyinference.CharacterStateMatrix;
import org.forester.phylogenyinference.CharacterStateMatrix.BinaryStates;

public class PaupLogParser {

    private static final String DATA_MATRIX_AND_RECONSTRUCTED_STATES_FOR_INTERNAL_NODES = "data matrix and reconstructed states for internal nodes";
    private Object              _nexus_source;

    private Object getNexusSource() {
        return _nexus_source;
    }

    public CharacterStateMatrix<BinaryStates> parse() throws IOException {
        final BufferedReader reader = ParserUtils.createReader( getNexusSource() );
        String line;
        boolean saw_line = false;
        int identifier_index = 0;
        boolean first_block = true;
        boolean saw_data_matrix_line = false;
        final List<String> identifiers = new ArrayList<String>();
        final List<List<BinaryStates>> states = new ArrayList<List<BinaryStates>>();
        boolean done = false;
        while ( ( ( line = reader.readLine() ) != null ) && !done ) {
            line = line.trim();
            if ( ( line.length() > 0 ) && !line.startsWith( "#" ) && !line.startsWith( ">" ) ) {
                if ( ( ( identifier_index > 0 ) && line.startsWith( "Tree " ) )
                        || line.startsWith( "Character change list" ) ) {
                    done = true;
                    continue;
                }
                if ( line.toLowerCase().startsWith( DATA_MATRIX_AND_RECONSTRUCTED_STATES_FOR_INTERNAL_NODES ) ) {
                    saw_line = false;
                    saw_data_matrix_line = true;
                    identifier_index = 0;
                    if ( first_block && ( line.indexOf( "continued" ) > 0 ) ) {
                        first_block = false;
                    }
                }
                if ( saw_data_matrix_line && line.startsWith( "----------" ) ) {
                    saw_line = true;
                }
                else if ( saw_line && ( line.indexOf( ' ' ) > 0 ) ) {
                    final String[] s = line.split( "\\s+" );
                    if ( s.length != 2 ) {
                        throw new NexusFormatException( "unexpected format at line: " + line );
                    }
                    final String identifier = s[ 0 ];
                    final String row = s[ 1 ];
                    if ( first_block ) {
                        if ( identifiers.contains( identifier ) ) {
                            throw new NexusFormatException( "identifier [" + identifier + "] is not unique in line: "
                                    + line );
                        }
                        identifiers.add( identifier );
                        states.add( new ArrayList<BinaryStates>() );
                    }
                    else {
                        if ( !identifiers.contains( identifier ) ) {
                            throw new NexusFormatException( "new identifier [" + identifier + "] at line: " + line );
                        }
                    }
                    for( int c = 0; c < row.length(); ++c ) {
                        final char ch = row.charAt( c );
                        if ( ch == '0' ) {
                            states.get( identifier_index ).add( BinaryStates.ABSENT );
                        }
                        else if ( ch == '1' ) {
                            states.get( identifier_index ).add( BinaryStates.PRESENT );
                        }
                        else {
                            throw new NexusFormatException( "unknown character state [" + ch + "] at line: " + line );
                        }
                    }
                    ++identifier_index;
                }
            }
        }
        final CharacterStateMatrix<BinaryStates> matrix = new BasicCharacterStateMatrix<BinaryStates>( states );
        int i = 0;
        for( final String identifier : identifiers ) {
            matrix.setIdentifier( i++, identifier );
        }
        return matrix;
    }

    public void setSource( final Object nexus_source ) throws PhylogenyParserException, IOException {
        if ( nexus_source == null ) {
            throw new PhylogenyParserException( getClass() + ": attempt to parse null object." );
        }
        _nexus_source = nexus_source;
    }
}
