# CourseNoteBeta

三门课程的讲义源码与 GitHub Pages 部署入口。

## 在线阅读

推送至 `main` 分支后，GitHub Actions 会渲染并发布站点：

**https://LZOI-ZJX.github.io/CourseNoteBeta/**

首页链接到三门课程子站点：`pde/`、`fuzzy/`、`topology/`。

## 仓库结构

| 目录 | 内容 |
|------|------|
| `pde/` | 偏微分方程数值解（Quarto Book：`lecture-notes/`、`Notes_MD/`） |
| `fuzzy/` | 模糊数学（Quarto Website，单页 `index.qmd`） |
| `topology/` | 拓扑学讲义翻译（Quarto Website） |

各子目录内含独立的 `_quarto.yml`，并已设置 GitHub Pages 子路径对应的 `website.site-url`。

## 本地预览

需安装 [Quarto](https://quarto.org/docs/get-started/)。

```bash
cd pde && quarto preview
```

```bash
cd fuzzy && quarto preview
```

```bash
cd topology && quarto preview
```

若需模拟线上子路径，可使用各 `_quarto.yml` 中的 `site-url` 进行校验。

## 许可

见仓库根目录 `LICENSE`（MIT）。
