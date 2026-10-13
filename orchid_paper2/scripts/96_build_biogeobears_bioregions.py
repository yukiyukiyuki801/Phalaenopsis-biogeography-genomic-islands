#!/usr/bin/env python3
"""Build biogeographic polygons and the 29-tip assignments used by BioGeoBEARS."""
from pathlib import Path
from osgeo import ogr, osr
import csv

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "19_biogeobears"
GADM = BASE / "inputs" / "gadm41"
OUT = BASE / "bioregions"
OUT.mkdir(parents=True, exist_ok=True)

SRS = osr.SpatialReference(); SRS.ImportFromEPSG(4326)

def read_geoms(code, level, predicate=lambda f: True):
    ds = ogr.Open(str(GADM / f"gadm41_{code}.gpkg"))
    lyr = ds.GetLayerByName(f"ADM_ADM_{level}")
    ans=[]
    for f in lyr:
        if predicate(f): ans.append(f.GetGeometryRef().Clone())
    return ans

def names(code, wanted):
    return read_geoms(code, 1, lambda f: f.GetField("NAME_1") in wanted)

def whole(code): return read_geoms(code, 0)

def union(gs):
    if not gs: raise ValueError("Cannot union empty geometry list")
    x=gs[0]
    for g in gs[1:]: x=x.Union(g)
    return x.MakeValid()

luzon={"Abra","Apayao","Aurora","Bataan","Batanes","Batangas","Benguet",
       "Bulacan","Cagayan","Cavite","Ifugao","Ilocos Norte","Ilocos Sur",
       "Isabela","Kalinga","La Union","Laguna","Metropolitan Manila",
       "Mountain Province","Nueva Ecija","Nueva Vizcaya","Pampanga",
       "Pangasinan","Quezon","Quirino","Rizal","Tarlac","Zambales"}
bicol={"Albay","Camarines Norte","Camarines Sur","Catanduanes","Masbate","Sorsogon"}
mindanao={"Agusan del Norte","Agusan del Sur","Basilan","Bukidnon","Camiguin",
          "Compostela Valley","Davao Oriental","Davao del Norte","Davao del Sur",
          "Dinagat Islands","Lanao del Norte","Lanao del Sur","Maguindanao",
          "Misamis Occidental","Misamis Oriental","North Cotabato","Sarangani",
          "South Cotabato","Sultan Kudarat","Sulu","Surigao del Norte",
          "Surigao del Sur","Tawi-Tawi","Zamboanga Sibugay","Zamboanga del Norte",
          "Zamboanga del Sur"}
greater_mindanao_islands={"Biliran","Bohol","Eastern Samar","Leyte",
                          "Northern Samar","Samar","Southern Leyte"}
greater_negros_panay={"Aklan","Antique","Capiz","Cebu","Guimaras","Iloilo",
                      "Masbate","Negros Occidental","Negros Oriental","Romblon",
                      "Siquijor"}
kalimantan={"Kalimantan Barat","Kalimantan Selatan","Kalimantan Tengah",
            "Kalimantan Timur","Kalimantan Utara"}
java={"Banten","Jakarta Raya","Jawa Barat","Jawa Tengah","Jawa Timur","Yogyakarta"}
sulawesi={"Bali","Gorontalo","Maluku Utara","Nusa Tenggara Barat",
          "Nusa Tenggara Timur","Sulawesi Barat","Sulawesi Selatan",
          "Sulawesi Tengah","Sulawesi Tenggara","Sulawesi Utara"}

# Operational biogeographic refinements: Kai (Maluku Tenggara/Tual) and the
# remaining Moluccas stay with Wallacea, whereas Aru is moved to the Sahul
# region with New Guinea. Queensland is clipped at 20 S to represent tropical
# northern Queensland rather than the entire state.
maluku_without_aru=read_geoms("IDN",2,lambda f: f.GetField("NAME_1")=="Maluku" and f.GetField("NAME_2")!="Kepulauan Aru")
aru=read_geoms("IDN",2,lambda f: f.GetField("NAME_1")=="Maluku" and f.GetField("NAME_2")=="Kepulauan Aru")
north_qld_clip=ogr.CreateGeometryFromWkt("POLYGON ((110 -20, 160 -20, 160 -8, 110 -8, 110 -20))")
qld_north=[g.Intersection(north_qld_clip).MakeValid() for g in names("AUS",{"Queensland"})]

# Luzon is divided by physiographic/endemism regions, not by a latitude line.
# Northern Luzon combines the Cordillera, Zambales and Northern Sierra Madre
# blocks plus the Batanes-Babuyan dispersal corridor. Southern Luzon combines
# the Southern Sierra Madre and Southern Luzon-Bicol volcanic belt. Province
# polygons are an operational approximation of these named mountain blocks:
# Aurora (including Casiguran) is retained with the Northern Sierra Madre,
# whereas Nueva Ecija (including Gabaldon) begins the southern/central side of
# the break. Zambales is retained north because it is a recognized northern
# centre of microendemism despite extending relatively far south.
luzon_complex=union(names("PHL",(luzon|bicol|{"Marinduque"})-{"Masbate"}))
north_luzon_provinces={"Abra","Apayao","Aurora","Batanes","Benguet",
    "Cagayan","Ifugao","Ilocos Norte","Ilocos Sur","Isabela","Kalinga",
    "La Union","Mountain Province","Nueva Vizcaya","Pangasinan","Quirino",
    "Zambales"}
north_luzon=union(names("PHL",north_luzon_provinces))
south_luzon=luzon_complex.Difference(north_luzon).MakeValid()

regions={
 "taiwan": union(whole("TWN")),
 "north_luzon": north_luzon,
 "south_luzon": south_luzon,
 "mindoro": union(names("PHL",{"Occidental Mindoro","Oriental Mindoro"})),
 "mindanao": union(names("PHL",mindanao|greater_mindanao_islands)),
 "greater_negros_panay": union(names("PHL",greater_negros_panay)),
 "palawan": union(names("PHL",{"Palawan"})),
 "borneo": union(names("IDN",kalimantan)+names("MYS",{"Sabah","Sarawak","Labuan"})+whole("BRN")),
 "java": union(names("IDN",java)),
 "sulawesi_moluccas_lesser_sunda": union(names("IDN",sulawesi)+maluku_without_aru+whole("TLS")),
 "new_guinea_queensland": union(names("IDN",{"Papua","Papua Barat"})+aru+whole("PNG")+qld_north),
}

gpkg=OUT/"biogeographic_regions.gpkg"
if gpkg.exists(): gpkg.unlink()
drv=ogr.GetDriverByName("GPKG"); ds=drv.CreateDataSource(str(gpkg))
lyr=ds.CreateLayer("bioregions",SRS,ogr.wkbMultiPolygon)
lyr.CreateField(ogr.FieldDefn("region_id",ogr.OFTString)); lyr.CreateField(ogr.FieldDefn("region_name",ogr.OFTString))
labels={
 "taiwan":"Taiwan","north_luzon":"Northern Luzon-Batanes-Babuyan",
 "south_luzon":"Southern Luzon-Bicol","mindoro":"Mindoro","mindanao":"Mindanao",
 "greater_negros_panay":"Greater Negros-Panay (unsampled)",
 "palawan":"Palawan","borneo":"Borneo","java":"Java",
 "sulawesi_moluccas_lesser_sunda":"Sulawesi-Moluccas-Lesser Sunda Islands",
 "new_guinea_queensland":"New Guinea-Queensland"}
for rid,g in regions.items():
    f=ogr.Feature(lyr.GetLayerDefn()); f.SetField("region_id",rid); f.SetField("region_name",labels[rid]); f.SetGeometry(g); lyr.CreateFeature(f)
ds=None

assign={
 "PAL":"palawan","PON":"borneo","SAB":"borneo","SK":"borneo","WJ":"java","SUR":"java",
 "SUL":"sulawesi_moluccas_lesser_sunda","AMB":"sulawesi_moluccas_lesser_sunda","TUA":"sulawesi_moluccas_lesser_sunda",
 "PNG":"new_guinea_queensland","QL":"new_guinea_queensland",
 "TW":"taiwan","FUG":"north_luzon","BAT":"north_luzon","ILO":"north_luzon",
 "APL":"north_luzon","AUR":"north_luzon","CAG":"north_luzon","APY":"north_luzon",
 "MAD":"north_luzon","CAS":"north_luzon","GAB":"south_luzon","BUE":"south_luzon",
 "MAR":"south_luzon","MDO":"mindoro","DO":"mindanao","LAM":"mindanao","MDA":"mindanao","PALp":"palawan"}
lineage={x:"PAM" for x in "PAL PON SAB SK WJ SUR SUL AMB TUA PNG QL".split()}
lineage.update({x:"PAP_PSA" for x in set(assign)-set(lineage)})
with open(OUT/"tip_to_thesis_bioregion.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["tip","lineage","region_id","region_name"])
    for tip in sorted(assign): w.writerow([tip,lineage[tip],assign[tip],labels[assign[tip]]])

print(f"Wrote {len(regions)} bioregions and {len(assign)} tip assignments to {OUT}")
