// $Id: BasicGoSubset.java,v 1.3 2009/11/10 19:57:09 cmzmasek Exp $
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

package org.forester.go;

public class BasicGoSubset implements GoSubset {

    final Type _type;

    public BasicGoSubset( final String s ) {
        final String my_s = s.trim().toLowerCase();
        if ( my_s.equals( GOSLIM_GENERIC_STR ) ) {
            _type = Type.GOSLIM_GENERIC;
        }
        else if ( my_s.equals( GOSLIM_GOA_STR ) ) {
            _type = Type.GOSLIM_GOA;
        }
        else if ( my_s.equals( GOSLIM_PIR_STR ) ) {
            _type = Type.GOSLIM_PIR;
        }
        else if ( my_s.equals( GOSUBSET_PROK_STR ) ) {
            _type = Type.GOSUBSET_PROK;
        }
        else if ( my_s.equals( GOSLIM_CANDIDA_STR ) ) {
            _type = Type.GOSLIM_CANDIDA;
        }
        else if ( my_s.equals( GOSLIM_PLANT_STR ) ) {
            _type = Type.GOSLIM_PLANT;
        }
        else if ( my_s.equals( GOSLIM_YEAST_STR ) ) {
            _type = Type.GOSLIM_YEAST;
        }
        else if ( my_s.equals( GOSLIM_POMBE_STR ) ) {
            _type = Type.GOSLIM_POMBE;
        }
        else {
            throw new IllegalArgumentException( "unknown GO subset type: " + my_s );
        }
    }

    public BasicGoSubset( final Type type ) {
        _type = type;
    }

    public int compareTo( final GoSubset sub ) {
        return getType().compareTo( sub.getType() );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check go subset equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check go subset equality to " + o + " [" + o.getClass()
                    + "]" );
        }
        else {
            return ( getType() == ( ( GoSubset ) o ).getType() );
        }
    }

    public Type getType() {
        return _type;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        switch ( getType() ) {
            case GOSLIM_CANDIDA:
                sb.append( GOSLIM_CANDIDA_STR );
                break;
            case GOSLIM_GENERIC:
                sb.append( GOSLIM_GENERIC_STR );
                break;
            case GOSLIM_GOA:
                sb.append( GOSLIM_GOA_STR );
                break;
            case GOSLIM_PIR:
                sb.append( GOSLIM_PIR_STR );
                break;
            case GOSLIM_PLANT:
                sb.append( GOSLIM_PLANT_STR );
                break;
            case GOSLIM_YEAST:
                sb.append( GOSLIM_YEAST_STR );
                break;
            case GOSUBSET_PROK:
                sb.append( GOSUBSET_PROK_STR );
                break;
            case GOSLIM_POMBE:
                sb.append( GOSLIM_POMBE_STR );
                break;
            default:
                new AssertionError( "unknown type: " + getType() );
        }
        return sb.toString();
    }
}
