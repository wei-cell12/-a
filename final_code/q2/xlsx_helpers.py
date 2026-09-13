"""问题二独立使用的轻量 XLSX 读取与模板样式工具。"""

from __future__ import annotations

import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


# 从指定 XLSX 工作表中读原始数据行。
def xlsx_rows(path: Path, sheet_index: int = 0) -> list[list[object]]:
    """按原始行读取 XLSX，不猜测表头。"""
    with zipfile.ZipFile(path) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", NS):
                shared.append("".join(node.text or "" for node in item.iter(f"{{{NS['m']}}}t")))
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        relation_map = {item.attrib["Id"]: item.attrib["Target"] for item in relationships}
        sheets = list(workbook.find("m:sheets", NS))
        relation_id = sheets[sheet_index].attrib[f"{{{NS['r']}}}id"]
        target = relation_map[relation_id].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        worksheet = ET.fromstring(archive.read(target))
        rows: list[list[object]] = []
        for row in worksheet.findall(".//m:sheetData/m:row", NS):
            values: dict[int, object] = {}
            for cell in row.findall("m:c", NS):
                reference = cell.attrib["r"]
                letters = "".join(ch for ch in reference if ch.isalpha())
                column = 0
                for letter in letters:
                    column = column * 26 + ord(letter.upper()) - 64
                value_node = cell.find("m:v", NS)
                cell_type = cell.attrib.get("t")
                value: object = None if value_node is None else value_node.text
                if cell_type == "s" and value is not None:
                    value = shared[int(value)]
                elif cell_type == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iter(f"{{{NS['m']}}}t"))
                elif value is not None:
                    value = float(value)
                values[column - 1] = value
            if values:
                rows.append([values.get(index) for index in range(max(values) + 1)])
        return rows


# 把从零开始的序号转换为 Excel 列字母。
def excel_column(index: int) -> str:
    letters = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


# 向 XLSX 样式表加入四位小数格式并返样式编号。
def add_four_decimal_style(styles_xml: bytes) -> tuple[bytes, int]:
    namespace = NS["m"]
    ET.register_namespace("", namespace)
    root = ET.fromstring(styles_xml)
    num_formats = root.find(f"{{{namespace}}}numFmts")
    if num_formats is None:
        num_formats = ET.Element(f"{{{namespace}}}numFmts", {"count": "1"})
        root.insert(0, num_formats)
    else:
        num_formats.attrib["count"] = str(int(num_formats.attrib.get("count", "0")) + 1)
    ET.SubElement(num_formats, f"{{{namespace}}}numFmt", {
        "numFmtId": "164", "formatCode": "0.0000",
    })
    cell_xfs = root.find(f"{{{namespace}}}cellXfs")
    if cell_xfs is None:
        raise ValueError("模板缺少 cellXfs 样式表")
    style_index = len(list(cell_xfs))
    style = ET.SubElement(cell_xfs, f"{{{namespace}}}xf", {
        "numFmtId": "164", "fontId": "0", "fillId": "0", "borderId": "0",
        "xfId": "0", "applyNumberFormat": "1", "applyAlignment": "1",
    })
    ET.SubElement(style, f"{{{namespace}}}alignment", {
        "horizontal": "center", "vertical": "center",
    })
    cell_xfs.attrib["count"] = str(style_index + 1)
    return ET.tostring(root, encoding="utf-8", xml_declaration=True), style_index
