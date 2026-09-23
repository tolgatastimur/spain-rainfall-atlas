# İspanya Yağış Atlası

İspanya'nın 1976'dan bugüne yağışını ERA5 verisiyle gösteren, hafif 3D bir web görselleştirmesi.

- **Harita:** Yarımada, özerk topluluk sınırları ve adlarıyla gösteriliyor. Balear (×1,6 büyütülmüş) ve Kanarya Adaları haritanın yanında ayrı kutularda duruyor. Yarımada 0,5°'lik (~50 km), adalar ise 0,25°'lik hücrelerle çiziliyor. Hücre değerleri Gauss yumuşatmasıyla tek parça bir yüzeye dönüşüyor. Tepelerin yüksekliği yağış miktarını, rengi seçilen ölçüyü (mm ya da 1991–2020 normaline oran) gösteriyor.
- **Zaman:** Slider ile istediğin döneme gidebilirsin. ▶ tuşu, slider'ın bulunduğu andan bugüne kadar **Gün / Hafta / Ay / Yıl** adımlarıyla oynatır (hız 1×, 4× veya 16×).
- **Yakınlaşma:** Bir topluluğa tıklayınca kamera oraya uçar, topluluğun sınırı çizilir ve yağış tepeleri yükselir. İller (Balear'da adalar) adlarıyla görünür, yan panelde o dönemin il bazında yağışı ve normale oranı listelenir. Bir ile tıklayınca aynı animasyonla il açılır. Üstteki yol (İspanya › Endülüs › Sevilla) veya Esc ile geri dönülür. Detay görünümleri yerel 0,25°'lik ERA5 noktalarını kullanır ve yalnızca açıldığında yüklenir.
- **Hikâye:** Veri bugüne ulaşınca (ya da "Su hikâyesi" düğmesine basınca) 6 sahnelik animasyonlu bölüm açılıyor. Sahneler şunlar: 50 yılın veriye dayalı özeti, su tüketimi, gelecekteki riskler, çözümler, kişisel tasarruf hesaplayıcısı, kaynaklar.

## Çalıştırma

```bash
python3 -m http.server 8765 --directory site
```

Ardından http://localhost:8765 adresini aç. Sayfa verileri `fetch` ile yüklediği için doğrudan dosya olarak (`file://`) açılmaz.

## Veri hattı

| Adım | Dosya | Açıklama |
|---|---|---|
| 1 | `pipeline/fetch_era5.py` | ERA5 saatlik toplam yağışı (0,25°) İspanya kutusu için indirir ve günlük toplamlara çevirir |
| 2 | `pipeline/build_geo.py` | Yarımadayı, özerk toplulukları ve adaları sadeleştirir |
| 3 | `pipeline/build_data.py` | 0,5° hücreler, 1991–2020 normali, yıllık toplamlar, `site/data/*` dosyaları |
| 4 | `pipeline/build_detail.py` | Her topluluk için 0,25°'lik noktalar, il/ada alt birimleri → `site/data/detail/*.json` |

```bash
python3 -m venv .venv && .venv/bin/pip install omfiles s3fs certifi numpy shapely
SSL_CERT_FILE=$(.venv/bin/python -m certifi) .venv/bin/python pipeline/fetch_era5.py 32
SSL_CERT_FILE=$(.venv/bin/python -m certifi) BOX=canarias .venv/bin/python pipeline/fetch_era5.py 32
.venv/bin/python pipeline/build_geo.py
.venv/bin/python pipeline/build_data.py
(cd pipeline && ../.venv/bin/python build_detail.py)
```

**Kaynak:** Copernicus Climate Change Service (C3S), ERA5 hourly data on single levels. Veri, Open-Meteo'nun açık S3 aynası üzerinden alınıyor (`s3://openmeteo/data/copernicus_era5/`). Bu ayna aynı ERA5 değerlerini API anahtarı gerektirmeden sunuyor. Veriyi güncellemek için 1. adımı (iki kutu için de) ve 3. adımı tekrar çalıştırman yeterli.

Sitedeki su kullanımı ve gelecek projeksiyonu rakamları (INE, MITECO, CEDEX, MedECC, IPCC) yuvarlanmış yaklaşık değerlerdir. Kaynakları hikâyenin son sahnesinde listelenmiştir.
