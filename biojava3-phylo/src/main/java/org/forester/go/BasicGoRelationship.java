// $Id: BasicGoRelationship.java,v 1.7 2009/01/13 19:49:29 cmzmasek Exp $
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

public class BasicGoRelationship implements GoRelationship {

    final Type _type;
    final GoId _go_id;

    public BasicGoRelationship( final String s ) {
        final String[] sa = s.split( " " );
        if ( sa.length != 2 ) {
            throw new IllegalArgumentException( "unexpected format for GO relationship: " + s );
        }
        final String type = sa[ 0 ].trim();
        final String go_id = sa[ 1 ].trim();
        if ( type.toLowerCase().equals( PART_OF_STR ) ) {
            _type = Type.PART_OF;
        }
        else if ( type.toLowerCase().equals( REGULATES_STR ) ) {
            _type = Type.REGULATES;
        }
        else if ( type.toLowerCase().equals( NEGATIVELY_REGULATES_STR ) ) {
            _type = Type.NEGATIVELY_REGULATES;
        }
        else if ( type.toLowerCase().equals( POSITIVELY_REGULATES_STR ) ) {
            _type = Type.POSITIVELY_REGULATES;
        }
        else {
            throw new IllegalArgumentException( "unknown GO relationship type: " + type );
        }
        _go_id = new GoId( go_id );
    }

    public BasicGoRelationship( final String type, final String go_id ) {
        if ( type.toLowerCase().equals( PART_OF_STR ) ) {
            _type = Type.PART_OF;
        }
        else {
            throw new IllegalArgumentException( "unknown GO relationship type: " + type );
        }
        _go_id = new GoId( go_id );
    }

    public BasicGoRelationship( final Type type, final GoId go_id ) {
        _type = type;
        _go_id = go_id;
    }

    public int compareTo( final GoRelationship rel ) {
        return getGoId().compareTo( rel.getGoId() );
    }

    /**
     * Based on value and type.
     * 
     * 
     */
    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check go relationship equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check go relationship equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return getType().equals( ( ( GoRelationship ) o ).getType() )
                    && getGoId().equals( ( ( GoRelationship ) o ).getGoId() );
        }
    }

    public GoId getGoId() {
        return _go_id;
    }

    public Type getType() {
        return _type;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        switch ( getType() ) {
            case PART_OF:
                sb.append( PART_OF_STR );
                break;
            case NEGATIVELY_REGULATES:
                sb.append( NEGATIVELY_REGULATES_STR );
                break;
            case POSITIVELY_REGULATES:
                sb.append( POSITIVELY_REGULATES_STR );
                break;
            case REGULATES:
                sb.append( REGULATES_STR );
                break;
            default:
                new AssertionError( "unknown type: " + getType() );
        }
        sb.append( ": " );
        sb.append( getGoId().toString() );
        return sb.toString();
    }
}
