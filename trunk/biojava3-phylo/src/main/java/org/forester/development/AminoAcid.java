// $Id: AminoAcid.java,v 1.6 2009/01/13 19:49:30 cmzmasek Exp $
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

package org.forester.development;

public final class AminoAcid {

    public static final char    UNKNOWN          = '?';
    public static final byte    UNKNOWN_CODE     = 20;
    public static final char    GAP              = '-';
    public static final byte    GAP_CODE         = 21;
    public static final char    TERMINATE        = '*';
    public static final byte    TERMINATE_CODE   = 22;
    private final static char[] AMINO_ACID_TABLE = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
            'F', 'P', 'S', 'T', 'W', 'Y', 'V', AminoAcid.UNKNOWN, AminoAcid.GAP, AminoAcid.TERMINATE };

    private AminoAcid() {
    }

    static public char getResidue( final byte state ) {
        return AminoAcid.AMINO_ACID_TABLE[ state ];
    }

    static public byte getState( final char c ) {
        switch ( c ) {
            case 'A':
            case 'a':
                return 0;
            case 'C':
            case 'c':
                return 4;
            case 'D':
            case 'd':
                return 3;
            case 'E':
            case 'e':
                return 6;
            case 'F':
            case 'f':
                return 13;
            case 'G':
            case 'g':
                return 7;
            case 'H':
            case 'h':
                return 8;
            case 'I':
            case 'i':
                return 9;
            case 'K':
            case 'k':
                return 11;
            case 'L':
            case 'l':
                return 10;
            case 'M':
            case 'm':
                return 12;
            case 'N':
            case 'n':
                return 2;
            case 'P':
            case 'p':
                return 14;
            case 'Q':
            case 'q':
                return 5;
            case 'R':
            case 'r':
                return 1;
            case 'S':
            case 's':
                return 15;
            case 'T':
            case 't':
                return 16;
            case 'V':
            case 'v':
                return 19;
            case 'W':
            case 'w':
                return 17;
            case 'Y':
            case 'y':
                return 18;
            case GAP:
                return AminoAcid.GAP_CODE;
            case TERMINATE:
                return AminoAcid.TERMINATE_CODE;
            default:
                return AminoAcid.UNKNOWN_CODE;
        }
    }

    static public boolean isGap( final byte state ) {
        return ( state == AminoAcid.GAP_CODE );
    }

    static public boolean isTerminate( final byte state ) {
        return ( state == AminoAcid.TERMINATE_CODE );
    }

    static public boolean isUnassignable( final byte state ) {
        return ( state > AminoAcid.TERMINATE_CODE );
    }

    static public boolean isUnknown( final byte state ) {
        return ( state == AminoAcid.UNKNOWN_CODE );
    }
}
