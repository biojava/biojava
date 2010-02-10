// $Id: Tuplet.java,v 1.2 2009/04/15 22:58:00 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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

package org.forester.sdi;

class Tuplet implements Comparable<Tuplet> {

    public static final int DEFAULT = -999;
    private final String    _key;
    private final double    _value1;
    private final double    _value2;
    private final double    _value3;
    private final double    _value4;
    private int[]           _p;            // Since

    Tuplet() {
        setSigns();
        _key = "";
        _value1 = Tuplet.DEFAULT;
        _value2 = Tuplet.DEFAULT;
        _value3 = Tuplet.DEFAULT;
        _value4 = Tuplet.DEFAULT;
    }

    // distance
    // needs to be
    // sorted in
    // different
    // direction than other values, and it is not
    // known which value will be the distance.
    Tuplet( final String name,
            final double value1,
            final double value2,
            final double value3,
            final double value4,
            final int c ) {
        setSigns();
        _key = name;
        _value1 = value1;
        _value2 = value2;
        _value3 = value3;
        _value4 = value4;
        if ( ( c >= 0 ) && ( c <= 3 ) ) {
            _p[ c ] = -1;
        }
    }

    Tuplet( final String name, final double value1, final double value2, final double value3, final int c ) {
        setSigns();
        _key = name;
        _value1 = value1;
        _value2 = value2;
        _value3 = value3;
        _value4 = Tuplet.DEFAULT;
        if ( ( c >= 0 ) && ( c <= 2 ) ) {
            _p[ c ] = -1;
        }
    }

    Tuplet( final String name, final double value1, final double value2, final int c ) {
        setSigns();
        _key = name;
        _value1 = value1;
        _value2 = value2;
        _value3 = Tuplet.DEFAULT;
        _value4 = Tuplet.DEFAULT;
        if ( ( c >= 0 ) && ( c <= 1 ) ) {
            _p[ c ] = -1;
        }
    }

    Tuplet( final String name, final double value1, final int c ) {
        setSigns();
        _key = name;
        _value1 = value1;
        _value2 = Tuplet.DEFAULT;
        _value3 = Tuplet.DEFAULT;
        _value4 = Tuplet.DEFAULT;
        if ( c == 0 ) {
            _p[ 0 ] = -1;
        }
    }

    public int compareTo( final Tuplet n ) {
        if ( ( getValue1() != Tuplet.DEFAULT ) && ( n.getValue1() != Tuplet.DEFAULT ) ) {
            if ( getValue1() < n.getValue1() ) {
                return _p[ 0 ];
            }
            if ( getValue1() > n.getValue1() ) {
                return ( -_p[ 0 ] );
            }
        }
        if ( ( getValue2() != Tuplet.DEFAULT ) && ( n.getValue2() != Tuplet.DEFAULT ) ) {
            if ( getValue2() < n.getValue2() ) {
                return _p[ 1 ];
            }
            if ( getValue2() > n.getValue2() ) {
                return ( -_p[ 1 ] );
            }
        }
        if ( ( getValue3() != Tuplet.DEFAULT ) && ( n.getValue3() != Tuplet.DEFAULT ) ) {
            if ( getValue3() < n.getValue3() ) {
                return _p[ 2 ];
            }
            if ( getValue3() > n.getValue3() ) {
                return ( -_p[ 2 ] );
            }
        }
        if ( ( getValue4() != Tuplet.DEFAULT ) && ( n.getValue4() != Tuplet.DEFAULT ) ) {
            if ( getValue4() < n.getValue4() ) {
                return _p[ 3 ];
            }
            if ( getValue4() > n.getValue4() ) {
                return ( -_p[ 3 ] );
            }
        }
        return ( getKey().compareTo( n.getKey() ) );
    }

    String getKey() {
        return _key;
    }

    double getValue1() {
        return _value1;
    }

    double getValue2() {
        return _value2;
    }

    double getValue3() {
        return _value3;
    }

    double getValue4() {
        return _value4;
    }

    private void setSigns() {
        _p = new int[ 4 ];
        _p[ 0 ] = _p[ 1 ] = _p[ 2 ] = _p[ 3 ] = +1;
    }
}