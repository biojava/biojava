// $Id: BasicGoXRef.java,v 1.11 2009/11/10 19:57:09 cmzmasek Exp $
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

public class BasicGoXRef implements GoXRef {

    private final String _xref;
    private final Type   _type;

    public BasicGoXRef( final String s ) {
        final String[] sa = s.split( ":" );
        if ( sa.length < 2 ) {
            throw new IllegalArgumentException( "unexpected format for GO xref: " + s );
        }
        final String type = sa[ 0 ].trim();
        if ( type.equals( EC_STR ) ) {
            _type = Type.EC;
        }
        else if ( type.equals( META_CYC_STR ) ) {
            _type = Type.META_CYC;
        }
        else if ( type.equals( REACTOME_STR ) ) {
            _type = Type.REACTOME;
        }
        else if ( type.equals( RESID_STR ) ) {
            _type = Type.RESID;
        }
        else if ( type.equals( UM_BBD_ENZYME_ID_STR ) ) {
            _type = Type.UM_BBD_ENZYME_ID;
        }
        else if ( type.equals( UM_BBD_PATHWAY_ID_STR ) ) {
            _type = Type.UM_BBD_PATHWAY_ID;
        }
        else if ( type.equals( UM_BBD_REACTIONID_STR ) ) {
            _type = Type.UM_BBD_REACTIONID;
        }
        else if ( type.equals( TC_STR ) ) {
            _type = Type.TC;
        }
        else if ( type.equals( ARACYC_STR ) ) {
            _type = Type.ARACYC;
        }
        else if ( type.equals( XX_STR ) ) {
            _type = Type.XX;
        }
        else if ( type.equals( PMID_STR ) ) {
            _type = Type.PMID;
        }
        else if ( type.equals( IMG_STR ) ) {
            _type = Type.IMG;
        }
        else if ( type.equals( GOC_STR ) ) {
            _type = Type.GOC;
        }
        else if ( type.equals( KEGG_STR ) ) {
            _type = Type.KEGG;
        }
        else if ( type.equals( WIKIPEDIA_STR ) ) {
            _type = Type.WIKIPEDIA;
        }
        else {
            throw new IllegalArgumentException( "unknown GO xref type: " + type );
        }
        _xref = sa[ 1 ].trim();
    }

    public BasicGoXRef( final Type type, final String xref ) {
        _type = type;
        _xref = xref;
    }

    public int compareTo( final GoXRef xref ) {
        return getXRef().compareTo( xref.getXRef() );
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
            throw new IllegalArgumentException( "attempt to check go xref equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check go xref equality to " + o + " [" + o.getClass() + "]" );
        }
        else {
            return getXRef().equals( ( ( GoXRef ) o ).getXRef() ) && getType().equals( ( ( GoXRef ) o ).getType() );
        }
    }

    public Type getType() {
        return _type;
    }

    @Override
    public String getXRef() {
        return _xref;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        switch ( getType() ) {
            case EC:
                sb.append( EC_STR );
                break;
            case META_CYC:
                sb.append( META_CYC_STR );
                break;
            case REACTOME:
                sb.append( REACTOME_STR );
                break;
            case RESID:
                sb.append( RESID_STR );
                break;
            case UM_BBD_ENZYME_ID:
                sb.append( UM_BBD_ENZYME_ID_STR );
                break;
            case UM_BBD_PATHWAY_ID:
                sb.append( UM_BBD_PATHWAY_ID_STR );
                break;
            case UM_BBD_REACTIONID:
                sb.append( UM_BBD_REACTIONID_STR );
                break;
            case TC:
                sb.append( TC_STR );
                break;
            case ARACYC:
                sb.append( ARACYC_STR );
                break;
            case XX:
                sb.append( XX_STR );
                break;
            case GOC:
                sb.append( GOC_STR );
                break;
            case IMG:
                sb.append( IMG_STR );
                break;
            case PMID:
                sb.append( PMID_STR );
                break;
            case WIKIPEDIA:
                sb.append( WIKIPEDIA_STR );
                break;
            default:
                new AssertionError( "unknown type: " + getType() );
        }
        sb.append( ":" );
        sb.append( getXRef() );
        return sb.toString();
    }
}
