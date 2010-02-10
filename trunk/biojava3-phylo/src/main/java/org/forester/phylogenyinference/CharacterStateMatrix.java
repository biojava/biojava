// $Id: CharacterStateMatrix.java,v 1.29 2009/02/04 01:52:18 cmzmasek Exp $
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

public interface CharacterStateMatrix<S> {

    public boolean containsCharacter( final String character );

    public boolean containsIdentifier( final String identifier );

    public CharacterStateMatrix<S> copy();

    public String getCharacter( int character_index );

    public int getCharacterIndex( final String character );

    public String getIdentifier( int identifier_index );

    public int getIdentifierIndex( final String identifier );

    public int getNumberOfCharacters();

    public int getNumberOfIdentifiers();

    public S getState( final int identifier_index, final int character_index );

    public S getState( final String identifier, final int character_index );

    public S getState( final String identifier, final String character );

    public boolean isEmpty();

    public CharacterStateMatrix<S> pivot();

    public void setCharacter( int character_index, final String character );

    public void setIdentifier( int identifier_index, final String identifier );

    public void setState( int identifier_index, int character_index, final S state );

    public void setState( final String identifier, int character_index, final S state );

    public void setState( final String identifier, final String character, final S state );

    public void toWriter( final Writer writer ) throws IOException;

    public void toWriter( final Writer writer, final Format format ) throws IOException;

    /**
     * It is crucial that the order
     * ABSENT, UNKNOWN, PRESENT not be changes since
     * this determines the sort order.
     *
     */
    static public enum BinaryStates {
        ABSENT, UNKNOWN, PRESENT;

        public char toChar() {
            switch ( this ) {
                case PRESENT:
                    return '1';
                case ABSENT:
                    return '0';
                case UNKNOWN:
                    return '?';
            }
            throw new IllegalStateException( "unknown state: " + this );
        }

        @Override
        public String toString() {
            switch ( this ) {
                case PRESENT:
                    return "1";
                case ABSENT:
                    return "0";
                case UNKNOWN:
                    return "?";
            }
            throw new IllegalStateException( "unknown state: " + this );
        }
    }

    public static enum Format {
        PHYLIP, FORESTER, NEXUS_BINARY
    }

    static public enum GainLossStates {
        GAIN, LOSS, UNCHANGED_PRESENT, UNCHANGED_ABSENT, UNKNOWN;

        @Override
        public String toString() {
            switch ( this ) {
                case GAIN:
                    return "+";
                case LOSS:
                    return "-";
                case UNCHANGED_PRESENT:
                    return "X";
                case UNCHANGED_ABSENT:
                    return ".";
                case UNKNOWN:
                    return "?";
            }
            throw new AssertionError( "unknown state: " + this );
        }
    }

    static public enum NucleotideStates {
        A, C, G, T, UNKNOWN;
    }
}
